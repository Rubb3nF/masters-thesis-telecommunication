% This script evaluates the steady-state tracking RMSE of two filters: the 
% Unscented Kalman Filter (UKF) and the WUPLKF.
% It performs a sweep over different SNR levels using DOA measurements 
% obtained from the motion_MUSIC estimator.

clear; close all; clc;


%% General & System Configuration 
% Filter & Time Parameters
kf_gen = pseudo_KF();
conversion_factor = 1000;          % Scale factor (km to meters)
sampling_period = 0.5;             % dt in seconds
num_samples = 150;                 % Total time steps
steady_state_start = 40;           % Step from which steady-state RMSE is calculated

% Process Noise Standard Deviations
accel_noise_std      = 0.1 * conversion_factor;  % m/s^2
accel_noise_std_PLKF = 0.01 * conversion_factor; % m/s^2

% SNR Sweep & Simulation Parameters
snr_range = -15:5:30;
snapshots = 500;
random_seed = 42;


%% Array Configuration 
num_elements_x = 16; 
num_elements_y = 16;
carrier_freq = 1e9;
element_spacing = 0.5;
light_speed = 3e8;

array_vca = UniformPlanarArray(num_elements_x, num_elements_y, element_spacing, element_spacing, carrier_freq, light_speed);
music_estimator = motion_MUSIC_opt(array_vca);

% Extract Sensor Positions
sensor_positions = array_vca.element_positions;
x_sensor = sensor_positions(1, :)';
y_sensor = sensor_positions(2, :)';

% MUSIC Search Grids
thetas_grid = linspace(-pi/2, 0, 91);
phis_grid   = linspace(-pi, pi, 361);


%% Scenario Definition (Trajectories) 
% Observer (Sensor) Trajectory
s_pos_hist = zeros(3, num_samples);
s_pos_hist(:,1) = [0; 60; 10] * conversion_factor; 
amplitude_lateral = 30 * conversion_factor;
angular_freq = 2 * pi / 40; 

for k = 2:num_samples
    t = (k-1) * sampling_period;
    pos_x = 1.0 * t * conversion_factor; 
    pos_y = s_pos_hist(2,1) + amplitude_lateral * sin(angular_freq * t);
    pos_z = s_pos_hist(3,1) + 3 * conversion_factor * sin(angular_freq * t / 2);
    s_pos_hist(:,k) = [pos_x; pos_y; pos_z];
end

% Target Ground Truth
target_real = zeros(6, num_samples);
target_real(1:3,1) = [80; 80; 1] * conversion_factor; 
target_real(4:6,1) = [-0.1; -0.2; 0.02] * conversion_factor;

% Constant Velocity State Transition Matrix
state_transition_matrix = [1, 0, 0, sampling_period, 0, 0;
                           0, 1, 0, 0, sampling_period, 0;
                           0, 0, 1, 0, 0, sampling_period; 
                           0, 0, 0, 1, 0, 0;
                           0, 0, 0, 0, 1, 0;
                           0, 0, 0, 0, 0, 1];

for k = 2:num_samples
    target_real(:,k) = state_transition_matrix * target_real(:,k-1);
end


%% Main SNR Sweep Loop
% Pre-allocate RMSE storage
num_snr = length(snr_range);
rmse_ukf_steady  = zeros(num_snr, 1);
rmse_plkf_steady = zeros(num_snr, 1);

for s_idx = 1:num_snr
    current_snr = snr_range(s_idx);
    fprintf('Processing SNR level: %2d dB (%d/%d)\n', current_snr, s_idx, num_snr);
    
    % Storage for current SNR measurements
    meas_u = zeros(3, num_samples);
    meas_angles = zeros(2, num_samples);
    sigma_angles_hist = zeros(2, num_samples);
    
    % Measurement Generation
    for k = 1:num_samples
        relative_pos = target_real(1:3, k) - s_pos_hist(:, k);
        dist_k = norm(relative_pos);
        true_phi = atan2(relative_pos(2), relative_pos(1));
        true_theta = asin(relative_pos(3) / dist_k); 
        
        % MUSIC DOA Estimation
        [theta_m, phi_m] = kf_gen.estimar_DOAs_MUSIC( ...
            num_elements_x, num_elements_y, element_spacing, carrier_freq, light_speed, ...
            true_theta, true_phi, snapshots, current_snr, false, ...
            thetas_grid, phis_grid, random_seed, music_estimator);
        
        % Coordinate Frame Adjustments
        theta_fixed = phi_m - pi;
        phi_fixed   = -theta_m;
        
        % Measurement Noise Covariance via CRB
        [var_th, var_ph] = crb_local(rad2deg(theta_m), rad2deg(phi_m), current_snr, snapshots, x_sensor, y_sensor);
        sigma_angles_hist(:, k) = [sqrt(var_th); sqrt(var_ph)];
        
        meas_angles(:, k) = [theta_fixed; phi_fixed];
        meas_u(:, k) = [cos(theta_fixed)*cos(phi_fixed); 
                        sin(theta_fixed)*cos(phi_fixed); 
                        sin(phi_fixed)];
    end
    
    % Filter Initialization 
    ukf_params.n = 6; ukf_params.alpha = 0.1; ukf_params.kappa = 0; ukf_params.beta = 2;
    
    q_pos = (sampling_period^4 / 4) * eye(3);
    q_pv  = (sampling_period^3 / 2) * eye(3);
    q_vel = (sampling_period^2)     * eye(3);
    
    process_q_PLKF = accel_noise_std_PLKF^2 * [q_pos, q_pv; q_pv, q_vel];
    process_q      = process_q_PLKF; % Shared process noise for comparison
    
    range_guess = 80 * conversion_factor;
    initial_state = [s_pos_hist(:,1) + range_guess * meas_u(:,1); 0; 0; 0];
    initial_cov   = diag([(100*conversion_factor)^2 * ones(1,3), (10*conversion_factor)^2 * ones(1,3)]);
    
    x_ukf  = initial_state; p_ukf  = initial_cov;
    x_plkf = initial_state; p_plkf = initial_cov;
    
    err_ukf_run  = zeros(1, num_samples);
    err_plkf_run = zeros(1, num_samples);
    
    % Tracking Loop 
    for k = 1:num_samples
        % UKF Processing
        z_ukf = meas_u(:, k);
        h_func = @(x) state_to_unitvector(x, s_pos_hist(:, k));
        
        th = meas_angles(1, k); 
        ph = meas_angles(2, k);
        J = [ -sin(th)*cos(ph), -cos(th)*sin(ph); 
               cos(th)*cos(ph), -sin(th)*sin(ph); 
               0,                cos(ph) ];
               
        r_ukf = J * diag(sigma_angles_hist(:, k)) * J';
        
        [x_ukf, p_ukf] = kf_gen.ukf_step(@(x) state_transition_matrix*x, h_func, x_ukf, p_ukf, z_ukf, process_q, r_ukf, ukf_params);
        err_ukf_run(k) = norm(x_ukf(1:3) - target_real(1:3, k));
        
        % PLKF Processing
        [x_plkf, p_plkf] = kf_gen.step_uplkf_3d(x_plkf, p_plkf, meas_angles(:, k), s_pos_hist(:, k), state_transition_matrix, process_q, sigma_angles_hist(:, k));
        err_plkf_run(k) = norm(x_plkf(1:3) - target_real(1:3, k));
    end
    
    % Steady-State Error Calculation 
    rmse_ukf_steady(s_idx)  = sqrt(mean(err_ukf_run(steady_state_start:end).^2));
    rmse_plkf_steady(s_idx) = sqrt(mean(err_plkf_run(steady_state_start:end).^2));
end


%% Visualization
% 3D Engagement Geometry
figure('Color', 'w', 'Name', 'Engagement Geometry', 'Position', [100 200 600 500]);
plot3(target_real(1,:), target_real(2,:), target_real(3,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Target Ground Truth'); hold on;
plot3(s_pos_hist(1,:), s_pos_hist(2,:), s_pos_hist(3,:), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Maneuvering Observer');
grid on; axis equal; view(45, 30);
xlabel('X Position (m)'); ylabel('Y Position (m)'); zlabel('Z Position (m)');
legend('Location', 'best');
title('3D Tracking Scenario Geometry');

% Performance Sensitivity (SNR vs RMSE)
figure('Color', 'w', 'Name', 'Performance Sensitivity', 'Position', [750 200 600 500]);
semilogy(snr_range, rmse_ukf_steady, 'ro-', 'LineWidth', 1.5, 'DisplayName', 'UKF RMSE'); hold on;
semilogy(snr_range, rmse_plkf_steady, 'bs--', 'LineWidth', 1.5, 'DisplayName', 'WUPLKF RMSE');
grid on;
xlabel('Signal to Noise Ratio (dB)'); 
ylabel('Steady State RMSE (m)');
legend('Location', 'best');
title('Estimation Accuracy Sensitivity to Measurement Noise');


%% LOCAL HELPER FUNCTIONS

function u = state_to_unitvector(x, s_pos)
    % Converts target state to a unit direction vector relative to the sensor
    diff_vec = x(1:3) - s_pos;
    magnitude = norm(diff_vec);
    if magnitude == 0
        u = [1; 0; 0]; 
    else
        u = diff_vec / magnitude; 
    end
end

function [var_theta, var_phi] = crb_local(theta_deg, phi_deg, snr_db, n_snapshots, px, py)
    % Computes the local Cramer-Rao Bound variances for a single source
    m_count = length(px); 
    theta_rad = deg2rad(theta_deg);
    phi_rad   = deg2rad(phi_deg);
    snr_linear = 10^(snr_db/10);
    
    u_coord = sin(theta_rad) * cos(phi_rad);
    v_coord = sin(theta_rad) * sin(phi_rad);
    steering_vec = exp(1j * 2 * pi * (px * u_coord + py * v_coord));
    
    du_dt = cos(theta_rad) * cos(phi_rad); 
    dv_dt = cos(theta_rad) * sin(phi_rad);
    du_dp = -sin(theta_rad) * sin(phi_rad); 
    dv_dp = sin(theta_rad) * cos(phi_rad);
    
    derivatives = [1j * 2 * pi * (px * du_dt + py * dv_dt) .* steering_vec, ...
                   1j * 2 * pi * (px * du_dp + py * dv_dp) .* steering_vec];
                   
    projection = eye(m_count) - (steering_vec * steering_vec') / m_count;
    fim_matrix = 2 * n_snapshots * snr_linear * real(derivatives' * projection * derivatives);
    
    crb_matrix = pinv(fim_matrix);
    var_theta = crb_matrix(1,1); 
    var_phi   = crb_matrix(2,2);
end