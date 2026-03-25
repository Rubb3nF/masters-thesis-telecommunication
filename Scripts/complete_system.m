% This script evaluates the tracking performance of the Unscented Kalman 
% Filter (UKF) against the Weighted Unbiased Pseudo-Linear Kalman Filter 
% (WUPLKF). It simulates a target performing a coordinated turn and an 
% observer moving in an S-curve to maintain high observability. Measurements 
% are generated using Motion-Aware MUSIC DOA estimation, with noise bounds derived 
% from the Cramer-Rao Bound.

clear; close all; clc;

%% General Configuration
kf_gen = pseudo_KF();
conversion_factor = 1000; % Convert km to meters
sampling_period = 0.5;    
num_samples = 190;        

accel_noise_std = 0.1 * conversion_factor;       % m/s^2 for UKF
accel_noise_std_PLKF = 0.01 * conversion_factor; % m/s^2 for WUPLKF

%% Array Configuration
num_elements_x = 16; 
num_elements_y = 16;
carrier_freq = 1e9;
element_spacing = 0.5;
light_speed = 3e8;

array_vca = UniformPlanarArray(num_elements_x, num_elements_y, ...
    element_spacing, element_spacing, carrier_freq, light_speed);
sensor_positions = array_vca.element_positions;
x_sensor = sensor_positions(1, :)';
y_sensor = sensor_positions(2, :)';

%% Observer Trajectory: S-Curve for High Observability 
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

%% Target Ground Truth: Constant Velocity
target_real = zeros(6, num_samples);
target_real(1:3,1) = [80; 80; 1] * conversion_factor;
target_real(4:6,1) = [-0.1; -0.1; 0.02] * conversion_factor;
% target_real(4:6,1) = [-0; -0; 0] * conversion_factor;
state_transition_matrix = [1 0 0 sampling_period 0 0;
0 1 0 0 sampling_period 0;
0 0 1 0 0 sampling_period;
0 0 0 1 0 0;
0 0 0 0 1 0;
0 0 0 0 0 1];
for k = 2:num_samples
target_real(:,k) = state_transition_matrix * target_real(:,k-1);
end

%% Target Ground Truth: Coordinated Turn
% target_real = zeros(6, num_samples);
% target_real(1:3,1) = [80; 80; 1] * conversion_factor; 
% 
% Initial velocity to start the movement
% v_init = 0.1 * conversion_factor; 
% target_real(4:6,1) = [-v_init; -v_init; 0.01 * conversion_factor]; 
% 
% Small turn rate (omega) to create a smooth curve
% omega = 0.02; 
% 
% for k = 2:num_samples
%     Rotation matrix for velocity in the XY plane
%     F_turn = [1, 0, 0,  sin(omega*sampling_period)/omega,      -(1-cos(omega*sampling_period))/omega, 0;
%               0, 1, 0, (1-cos(omega*sampling_period))/omega,  sin(omega*sampling_period)/omega,      0;
%               0, 0, 1,  0,                                     0,                                     sampling_period;
%               0, 0, 0,  cos(omega*sampling_period),           -sin(omega*sampling_period),            0;
%               0, 0, 0,  sin(omega*sampling_period),            cos(omega*sampling_period),            0;
%               0, 0, 0,  0,                                     0,                                     1];
%     target_real(:,k) = F_turn * target_real(:,k-1);
% end

%% Measurement Generation: MUSIC and CRB 
music_estimator = motion_MUSIC_opt(array_vca);
thetas_grid = linspace(-pi/2, 0, 91);
phis_grid   = linspace(-pi, pi, 361);
snapshots = 500;
snr_val = 15;
random_seed = 42;

meas_u = zeros(3, num_samples);
meas_angles = zeros(2, num_samples);
sigma_angles_hist = zeros(2, num_samples);

for k = 1:num_samples
    relative_pos = target_real(1:3, k) - s_pos_hist(:, k);
    dist_k = norm(relative_pos);
    
    true_phi = atan2(relative_pos(2), relative_pos(1));
    true_theta = asin(relative_pos(3) / dist_k); 
    
    [theta_m, phi_m] = kf_gen.estimar_DOAs_MUSIC( ...
        num_elements_x, num_elements_y, element_spacing, carrier_freq, light_speed, ...
        true_theta, true_phi, snapshots, snr_val, false, ...
        thetas_grid, phis_grid, random_seed, music_estimator);
    
    % Adjust estimated angles to tracking coordinate system
    theta_fixed = phi_m - pi;
    phi_fixed = -theta_m;
    
    % Compute Cramer-Rao Bound for noise variance
    [v_th, v_ph] = crb_local(rad2deg(phi_m), rad2deg(theta_m), snr_val, snapshots, x_sensor, y_sensor);
    sigma_angles_hist(:, k) = [sqrt(v_th); sqrt(v_ph)];
    
    meas_angles(:, k) = [theta_fixed; phi_fixed];
    meas_u(:, k) = [cos(theta_fixed)*cos(phi_fixed); sin(theta_fixed)*cos(phi_fixed); sin(phi_fixed)];
end

%% Shared Filter Initialization
ukf_params.n = 6;
ukf_params.alpha = 0.1;
ukf_params.kappa = 0;
ukf_params.beta  = 2;

q_pos = (sampling_period^4 / 4) * eye(3);
q_pv  = (sampling_period^3 / 2) * eye(3);
q_vel = (sampling_period^2) * eye(3);

process_q_PLKF = accel_noise_std_PLKF^2 * [q_pos, q_pv; q_pv, q_vel];
process_q = process_q_PLKF; % Shared process noise for comparison

range_guess = 80 * conversion_factor;
initial_state = [s_pos_hist(:,1) + range_guess * meas_u(:,1); 0; 0; 0];
initial_cov = diag([(100*conversion_factor)^2*ones(1,3), (10*conversion_factor)^2*ones(1,3)]);

% State holders
x_ukf  = initial_state; p_ukf  = initial_cov;
x_plkf = initial_state; p_plkf = initial_cov;

ukf_est_hist  = zeros(6, num_samples);
plkf_est_hist = zeros(6, num_samples);

%% Unified Estimation Loop 
state_transition_matrix = [1, 0, 0, sampling_period, 0, 0;
                           0, 1, 0, 0, sampling_period, 0;
                           0, 0, 1, 0, 0, sampling_period; 
                           0, 0, 0, 1, 0, 0;
                           0, 0, 0, 0, 1, 0;
                           0, 0, 0, 0, 0, 1];

for k = 1:num_samples
    
    % UKF Filter Step
    z_ukf = meas_u(:, k);
    h_func = @(x) state_to_unitvector(x, s_pos_hist(:, k));
    
    th = meas_angles(1, k); 
    ph = meas_angles(2, k);
    
    % Jacobian to transform angle noise to unit vector space
    J = [ -sin(th)*cos(ph), -cos(th)*sin(ph);
           cos(th)*cos(ph), -sin(th)*sin(ph);
           0,                cos(ph) ];
                
    r_ukf = J * diag(sigma_angles_hist(:, k)) * J';
    [x_ukf, p_ukf] = kf_gen.ukf_step(@(x) state_transition_matrix*x, h_func, ...
                                     x_ukf, p_ukf, z_ukf, process_q, r_ukf, ukf_params);
    ukf_est_hist(:, k) = x_ukf;
   
    % WUPLKF Filter Step
    [x_plkf, p_plkf] = kf_gen.step_uplkf_3d(x_plkf, p_plkf, ...
                                            meas_angles(:, k), ...
                                            s_pos_hist(:, k), ...
                                            state_transition_matrix, ...
                                            process_q_PLKF, ...
                                            sigma_angles_hist(:, k));
    plkf_est_hist(:, k) = x_plkf;
end

%% Visualization
% Trajectories 
figure('Color','w','Name','Tracking Comparison', 'Position', [100, 100, 800, 600]);
plot3(target_real(1,:), target_real(2,:), target_real(3,:), 'k-', 'LineWidth', 2.5); hold on;
plot3(ukf_est_hist(1,:), ukf_est_hist(2,:), ukf_est_hist(3,:), 'r--', 'LineWidth', 1.5);
plot3(plkf_est_hist(1,:), plkf_est_hist(2,:), plkf_est_hist(3,:), 'm:', 'LineWidth', 2);
plot3(s_pos_hist(1,:), s_pos_hist(2,:), s_pos_hist(3,:), 'b:', 'LineWidth', 1.5);

grid on; axis equal; view(45, 30);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend('Ground Truth', 'UKF', 'WUPLKF', 'Observer', 'Location', 'best');
title('Integrated 3D Target Tracking Performance');

% Euclidean Error 
err_ukf  = sqrt(sum((ukf_est_hist(1:3,:)  - target_real(1:3,:)).^2, 1));
err_plkf = sqrt(sum((plkf_est_hist(1:3,:) - target_real(1:3,:)).^2, 1));

figure('Color','w','Name','Error Comparison', 'Position', [950, 100, 600, 400]);
plot(err_ukf, 'r', 'LineWidth', 1.5); hold on;
plot(err_plkf, 'm', 'LineWidth', 1.5);
grid on; 
xlabel('Time Step'); ylabel('Euclidean Error (m)');
legend('UKF Error', 'WUPLKF Error', 'Location', 'best');
title('Tracking Error Performance Comparison');

% Calculate mean steady-state error (from step 85 onwards)
steady_state_idx = 85;
mean_err_ukf  = mean(err_ukf(steady_state_idx:end));
mean_err_plkf = mean(err_plkf(steady_state_idx:end));

fprintf('\n--- Steady-State Tracking Error (Steps %d to %d) ---\n', steady_state_idx, num_samples);
fprintf('Average Position Error (UKF):    %.4f m\n', mean_err_ukf);
fprintf('Average Position Error (WUPLKF): %.4f m\n', mean_err_plkf);


%% LOCAL HELPER FUNCTIONS

function u = state_to_unitvector(x, s_pos)
    % Converts an absolute Cartesian state vector to an observer-relative unit vector
    diff_vec = x(1:3) - s_pos;
    magnitude = norm(diff_vec);
    if magnitude == 0
        u = [1; 0; 0];
    else
        u = diff_vec / magnitude;
    end
end

function [var_theta, var_phi] = crb_local(theta_deg, phi_deg, snr_db, n_snapshots, px, py)
    % Computes the local Cramer-Rao Bound (CRB) for a single source DOA estimation
    m_count = length(px); 
    theta_rad = deg2rad(theta_deg);
    phi_rad   = deg2rad(phi_deg);
    snr_linear = 10^(snr_db/10);
    
    u_coord = sin(theta_rad)*cos(phi_rad);
    v_coord = sin(theta_rad)*sin(phi_rad);
    steering_vec = exp(1j * 2 * pi * (px * u_coord + py * v_coord));
    
    du_dt = cos(theta_rad)*cos(phi_rad);
    dv_dt = cos(theta_rad)*sin(phi_rad);
    du_dp = -sin(theta_rad)*sin(phi_rad);
    dv_dp =  sin(theta_rad)*cos(phi_rad);
    
    derivatives = [1j*2*pi*(px*du_dt + py*dv_dt).*steering_vec, ...
                   1j*2*pi*(px*du_dp + py*dv_dp).*steering_vec];
         
    projection = eye(m_count) - (steering_vec * steering_vec') / m_count;
    fim_matrix = 2 * n_snapshots * snr_linear * real(derivatives' * projection * derivatives);
    crb_matrix = pinv(fim_matrix);
    
    var_theta = crb_matrix(1,1);
    var_phi   = crb_matrix(2,2);
end