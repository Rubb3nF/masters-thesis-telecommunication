%% Unified Multi Target Tracking Comparison UKF vs WUPLKF
% This script evaluates the steady state tracking RMSE of two filters the 
% Unscented Kalman Filter UKF and the WUPLKF
% It performs a sweep over different SNR levels using DOA measurements 
% obtained from the motion MUSIC estimator

clear; close all; clc;

%% General Configuration
kf_gen = pseudo_KF();
conversion_factor = 1000; 
sampling_period = 0.5;    
num_samples = 190;        
accel_noise_std = 0.1 * conversion_factor; % m/s^2
accel_noise_std_PLKF = 0.01 * conversion_factor; % m/s^2

%% Array Configuration
num_elements_x = 16; 
num_elements_y = 16;
carrier_freq = 1e9;
element_spacing = 0.5;
light_speed = 3e8;
array_vca = UniformPlanarArray(num_elements_x, num_elements_y, element_spacing, element_spacing, carrier_freq, light_speed);
sensor_positions = array_vca.element_positions;
x_sensor = sensor_positions(1, :)';
y_sensor = sensor_positions(2, :)';

%% Observer Trajectory S Curve for High Observability
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

%% Target Ground Truth Constant Velocity 2 Targets
% Target Right to left
target1_real = zeros(6, num_samples);
target1_real(1:3,1) = [80; 80; 1] * conversion_factor; 
target1_real(4:6,1) = [-0.1; -0.1; 0.02] * conversion_factor; 

% Target Left to right
target2_real = zeros(6, num_samples);
target2_real(1:3,1) = [20; 80; 1] * conversion_factor; 
target2_real(4:6,1) = [0.15; -0.1; 0.01] * conversion_factor; 

state_transition_matrix = [1 0 0 sampling_period 0 0;
                           0 1 0 0 sampling_period 0;
                           0 0 1 0 0 sampling_period; 
                           0 0 0 1 0 0;
                           0 0 0 0 1 0;
                           0 0 0 0 0 1];
for k = 2:num_samples
    target1_real(:,k) = state_transition_matrix * target1_real(:,k-1);
    target2_real(:,k) = state_transition_matrix * target2_real(:,k-1);
end

%% Measurement Generation MUSIC and CRB Multi Target
signal_model = SignalModel(array_vca);
music_estimator = motion_MUSIC_opt(array_vca);
snapshots = 500;
snr_val = 15;

% Memory preallocation for measurements
meas_u1 = zeros(3, num_samples); meas_angles1 = zeros(2, num_samples); sigma_angles_hist1 = zeros(2, num_samples);
meas_u2 = zeros(3, num_samples); meas_angles2 = zeros(2, num_samples); sigma_angles_hist2 = zeros(2, num_samples);
last_meas = zeros(2, 2); 

for k = 1:num_samples
    % True DOAs
    rel1 = target1_real(1:3, k) - s_pos_hist(:, k);
    true_ph1 = atan2(rel1(2), rel1(1));
    true_th1 = asin(rel1(3) / norm(rel1)); 
    
    rel2 = target2_real(1:3, k) - s_pos_hist(:, k);
    true_ph2 = atan2(rel2(2), rel2(1));
    true_th2 = asin(rel2(3) / norm(rel2));
    
    % Signal Generation and Estimation
    [X, ~, ~] = signal_model.generateSignals([true_th1, true_th2], [true_ph1, true_ph2], ...
                snapshots, snr_val, 'complex_sinusoid', 0, 'white', k);
    [th_raw, ph_raw] = music_estimator.estimate(X, 2, 0.1*pi/180, true);
    
    % Simple Hungarian Association
    raw_m = [th_raw(:).'; ph_raw(:).'];
    if k == 1
        [~, idx] = sort(raw_m(2,:)); % Sort by initial azimuth
        meas = raw_m(:, idx);
    else
        dists = pdist2(raw_m', last_meas');
        if dists(1,1) + dists(2,2) < dists(1,2) + dists(2,1)
            meas = raw_m;
        else
            meas = raw_m(:, [2 1]);
        end
    end
    last_meas = meas;
    
    % TARGET Processing Measurement
    th_m1 = meas(1,1); ph_m1 = meas(2,1);
    th_fixed1 = ph_m1 - pi; ph_fixed1 = -th_m1;
    [v_th1, v_ph1] = crb_local(rad2deg(ph_m1), rad2deg(th_m1), snr_val, snapshots, x_sensor, y_sensor);
    
    meas_angles1(:, k) = [th_fixed1; ph_fixed1];
    sigma_angles_hist1(:, k) = [sqrt(v_th1); sqrt(v_ph1)];
    meas_u1(:, k) = [cos(th_fixed1)*cos(ph_fixed1); sin(th_fixed1)*cos(ph_fixed1); sin(ph_fixed1)];
    
    % TARGET Processing Measurement
    th_m2 = meas(1,2); ph_m2 = meas(2,2);
    th_fixed2 = ph_m2 - pi; ph_fixed2 = -th_m2;
    [v_th2, v_ph2] = crb_local(rad2deg(ph_m2), rad2deg(th_m2), snr_val, snapshots, x_sensor, y_sensor);
    
    meas_angles2(:, k) = [th_fixed2; ph_fixed2];
    sigma_angles_hist2(:, k) = [sqrt(v_th2); sqrt(v_ph2)];
    meas_u2(:, k) = [cos(th_fixed2)*cos(ph_fixed2); sin(th_fixed2)*cos(ph_fixed2); sin(ph_fixed2)];
end

%% Shared Filter Initialization
ukf_params.n = 6;
ukf_params.alpha = 0.1;
ukf_params.kappa = 0;
ukf_params.beta  = 2;
q_pos = (sampling_period^4 / 4) * eye(3);
q_pv  = (sampling_period^3 / 2) * eye(3);
q_vel = (sampling_period^2) * eye(3);
process_q = accel_noise_std^2 * [q_pos, q_pv; q_pv, q_vel];
process_q_PLKF = accel_noise_std_PLKF^2 * [q_pos, q_pv; q_pv, q_vel];
range_guess = 80 * conversion_factor;
initial_cov = diag([(100*conversion_factor)^2*ones(1,3), (10*conversion_factor)^2*ones(1,3)]);

% Target 1 Initialization
x_ukf1 = [s_pos_hist(:,1) + range_guess * meas_u1(:,1); 0; 0; 0]; p_ukf1 = initial_cov;
x_plkf1 = x_ukf1; p_plkf1 = initial_cov;
ukf_est_hist1 = zeros(6, num_samples); plkf_est_hist1 = zeros(6, num_samples);

% Target 2 Initialization
x_ukf2 = [s_pos_hist(:,1) + range_guess * meas_u2(:,1); 0; 0; 0]; p_ukf2 = initial_cov;
x_plkf2 = x_ukf2; p_plkf2 = initial_cov;
ukf_est_hist2 = zeros(6, num_samples); plkf_est_hist2 = zeros(6, num_samples);

%% Unified Estimation Loop
for k = 1:num_samples
    h_func = @(x) state_to_unitvector(x, s_pos_hist(:, k));
    
    % TARGET 1
    % UKF 1
    th1 = meas_angles1(1, k); ph1 = meas_angles1(2, k);
    J1 = [ -sin(th1)*cos(ph1), -cos(th1)*sin(ph1);
            cos(th1)*cos(ph1), -sin(th1)*sin(ph1);
            0,                  cos(ph1) ];
    r_ukf1 = J1 * diag(sigma_angles_hist1(:, k)) * J1';
    
    [x_ukf1, p_ukf1] = kf_gen.ukf_step(@(x) state_transition_matrix*x, h_func, x_ukf1, p_ukf1, meas_u1(:, k), process_q, r_ukf1, ukf_params);
    ukf_est_hist1(:, k) = x_ukf1;
   
    % PLKF 1 WUPLKF
    [x_plkf1, p_plkf1] = kf_gen.step_uplkf_3d(x_plkf1, p_plkf1, meas_angles1(:, k), s_pos_hist(:, k), ...
                                              state_transition_matrix, process_q_PLKF, sigma_angles_hist1(:, k));
    plkf_est_hist1(:, k) = x_plkf1;
    
    % TARGET 2
    % UKF 2
    th2 = meas_angles2(1, k); ph2 = meas_angles2(2, k);
    J2 = [ -sin(th2)*cos(ph2), -cos(th2)*sin(ph2);
            cos(th2)*cos(ph2), -sin(th2)*sin(ph2);
            0,                  cos(ph2) ];
    r_ukf2 = J2 * diag(sigma_angles_hist2(:, k)) * J2';
    
    [x_ukf2, p_ukf2] = kf_gen.ukf_step(@(x) state_transition_matrix*x, h_func, x_ukf2, p_ukf2, meas_u2(:, k), process_q, r_ukf2, ukf_params);
    ukf_est_hist2(:, k) = x_ukf2;
   
    % PLKF 2 WUPLKF
    [x_plkf2, p_plkf2] = kf_gen.step_uplkf_3d(x_plkf2, p_plkf2, meas_angles2(:, k), s_pos_hist(:, k), ...
                                              state_transition_matrix, process_q_PLKF, sigma_angles_hist2(:, k));
    plkf_est_hist2(:, k) = x_plkf2;
end

%% Visualization Trajectories
figure('Color','w','Name','Tracking Comparison', 'Position', [100 100 800 600]);
hold on; grid on; axis equal; view(45, 30);
% Observer
plot3(s_pos_hist(1,:), s_pos_hist(2,:), s_pos_hist(3,:), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Observer');
% Target 1
plot3(target1_real(1,:), target1_real(2,:), target1_real(3,:), 'k-', 'LineWidth', 2.5, 'DisplayName', 'T1 Ground Truth');
plot3(ukf_est_hist1(1,:), ukf_est_hist1(2,:), ukf_est_hist1(3,:), 'g--', 'LineWidth', 1.5, 'DisplayName', 'T1 UKF');
plot3(plkf_est_hist1(1,:), plkf_est_hist1(2,:), plkf_est_hist1(3,:), 'm:', 'LineWidth', 1.5, 'DisplayName', 'T1 WUPLKF');
% Target 2
plot3(target2_real(1,:), target2_real(2,:), target2_real(3,:), 'k-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
plot3(ukf_est_hist2(1,:), ukf_est_hist2(2,:), ukf_est_hist2(3,:), 'b--', 'LineWidth', 1.5, 'DisplayName', 'T2 UKF');
plot3(plkf_est_hist2(1,:), plkf_est_hist2(2,:), plkf_est_hist2(3,:), 'c:', 'LineWidth', 1.5, 'DisplayName', 'T2 WUPLKF');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend('Location', 'best');
title('Integrated 3D Multi Target Tracking');

%% Visualization Euclidean Error
err_ukf1 = sqrt(sum((ukf_est_hist1(1:3,:) - target1_real(1:3,:)).^2, 1));
err_plkf1 = sqrt(sum((plkf_est_hist1(1:3,:) - target1_real(1:3,:)).^2, 1));
err_ukf2 = sqrt(sum((ukf_est_hist2(1:3,:) - target2_real(1:3,:)).^2, 1));
err_plkf2 = sqrt(sum((plkf_est_hist2(1:3,:) - target2_real(1:3,:)).^2, 1));

figure('Color','w','Name','Error Comparison', 'Position', [950 100 600 600]);
subplot(2,1,1);
semilogy(err_ukf1, 'g', 'LineWidth', 1.5); hold on;
semilogy(err_plkf1, 'm', 'LineWidth', 1.5);
grid on; ylabel('Error (m)'); title('Target 1 Error'); legend('UKF', 'WUPLKF');
subplot(2,1,2);
semilogy(err_ukf2, 'b', 'LineWidth', 1.5); hold on;
semilogy(err_plkf2, 'c', 'LineWidth', 1.5);
grid on; xlabel('Time Step'); ylabel('Error (m)'); title('Target 2 Error'); legend('UKF', 'WUPLKF');

% Steady State Calculations Assuming sample 85 to the end
ss_start = 85;
fprintf('\n--- Average Position Error (Samples %d to %d) ---\n', ss_start, num_samples);
fprintf('Target 1 -> UKF: %.4f m | WUPLKF: %.4f m\n', mean(err_ukf1(ss_start:end)), mean(err_plkf1(ss_start:end)));
fprintf('Target 2 -> UKF: %.4f m | WUPLKF: %.4f m\n', mean(err_ukf2(ss_start:end)), mean(err_plkf2(ss_start:end)));

%% Helper Functions
function u = state_to_unitvector(x, s_pos)
    diff_vec = x(1:3) - s_pos;
    magnitude = norm(diff_vec);
    if magnitude == 0
        u = [1; 0; 0];
    else
        u = diff_vec / magnitude;
    end
end

function [var_theta, var_phi] = crb_local(theta_deg, phi_deg, snr_db, n_snapshots, px, py)
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