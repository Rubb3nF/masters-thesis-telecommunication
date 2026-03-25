% This script runs a Monte Carlo simulation over a sweep of Signal to Noise 
% Ratios to validate the DOA estimation performance of Motion-Aware MUSIC 
% against the theoretical Cramer-Rao Bound. 
% 
% It evaluates the percentage of estimates falling within the 
% confidence ellipse, derived from the Fisher Information Matrix:
% FIM = 2 * N_snap * SNR * real(D' * P_perp * D)

clear; clc; close all;

%% System Parameters
% Array & Signal Parameters
M = 16; 
N = 16; 
d = 0.5; 
fc = 1e9; 
c = 3e8; 

% Spatial Smoothing Sub-array Configuration
M_sub = 16; 
N_sub = 16; 

% Initialize Objects
array_vca = UniformPlanarArray(M, N, d, d, fc, c);
signal_model = SignalModel(array_vca);
estimator = motion_MUSIC_opt(array_vca);

% Sensor Positions
positions = array_vca.element_positions;
x_pos = positions(1, :)';
y_pos = positions(2, :)';

% Target True Position
theta_true_deg = 69.99;
phi_true_deg   = -14.04;
theta_true = deg2rad(theta_true_deg);
phi_true   = deg2rad(phi_true_deg);

% Simulation Setup
N_snapshots = 300;
num_montecarlo = 500;
sigma_factor = 1; 

% SNR Sweep Setup
snr_db_sweep = -20:4:20;
num_snr = length(snr_db_sweep);
percentage_inside_list = zeros(num_snr, 1);


%% Sub-array Geometry Initialization
vec_x_sub = (0:M_sub-1) - (M_sub-1)/2;
vec_y_sub = (0:N_sub-1) - (N_sub-1)/2;
[X_sub, Y_sub] = ndgrid(vec_x_sub, vec_y_sub);
p_x_sub = X_sub(:) * d;
p_y_sub = Y_sub(:) * d;


%% Monte Carlo Loop Over SNR
wb = waitbar(0, 'Running SNR sweep');

for snr_idx = 1:num_snr
    snr_db = snr_db_sweep(snr_idx);
    
    errors_theta = zeros(num_montecarlo, 1);
    errors_phi   = zeros(num_montecarlo, 1);
    
    for mc = 1:num_montecarlo
        % Update waitbar smoothly
        waitbar((snr_idx - 1 + mc/num_montecarlo) / num_snr, wb, ...
                sprintf('SNR: %d dB (Iteration %d/%d)', snr_db, mc, num_montecarlo));
        
        seed = 100 + mc;
        
        % Generate Signal
        [X, ~, ~] = signal_model.generateSignals( ...
            theta_true, phi_true, ...
            N_snapshots, snr_db, ...
            'complex_sinusoid', 0.0, 'white', seed);
        
        % Initial Guess (Grid Search)
        theta_grid = linspace(theta_true - 0.1, theta_true + 0.1, 20);
        phi_grid   = linspace(phi_true - 0.1, phi_true + 0.1, 20);
        
        try
            [theta_init, phi_init] = estimator.estimate(X, 1, theta_grid, phi_grid, true);
            start_point = [theta_init, phi_init];
        catch
            start_point = [theta_true, phi_true];
        end
        
        % Subspace Processing & Spatial Smoothing
        R_ss = spatial_smoothing(X, M, N, M_sub, N_sub);
        [E, Dvals] = eig(R_ss);
        [~, idx] = sort(diag(Dvals), 'ascend');
        En = E(:, idx(1:end-1));
        UU = En * En';
        
        % High-Resolution Refinement
        cost_func = @(ang) music_cost_function(ang, p_x_sub, p_y_sub, UU);
        options = optimset('Display', 'off', 'TolX', 1e-9);
        refined_pos = fminsearch(cost_func, start_point, options);
        
        errors_theta(mc) = refined_pos(1) - theta_true;
        errors_phi(mc)   = refined_pos(2) - phi_true;
    end
    
    % Compute Theoretical CRB for current SNR
    [var_theta_crb, var_phi_crb] = crb( ...
        theta_true_deg, phi_true_deg, snr_db, ...
        N_snapshots, x_pos, y_pos);
    
    % Check how many estimates fall within the confidence ellipse
    inside = (errors_theta.^2 ./ (sigma_factor^2 * var_theta_crb) + ...
              errors_phi.^2   ./ (sigma_factor^2 * var_phi_crb)) <= 1;
              
    percentage_inside_list(snr_idx) = 100 * sum(inside) / num_montecarlo;
end
close(wb);


%% Visualization
% SNR vs Percentage Inside Region
figure('Color', 'w', 'Name', 'Confidence Region Validation', 'Position', [100 200 600 450]);
plot(snr_db_sweep, percentage_inside_list, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
xlabel('SNR (dB)');
ylabel(sprintf('Percentage Inside %d-\\sigma CRB Region', sigma_factor));
title(sprintf('CRB %d-\\sigma Confidence Region Validation vs SNR', sigma_factor));
grid on; ylim([0 105]);


% Plot 2: Error Distribution for Last SNR 
figure('Color', 'w', 'Name', 'Error Scatter Plot', 'Position', [750 200 600 500]);
plot(rad2deg(errors_theta), rad2deg(errors_phi), 'bx', 'LineWidth', 1.2, 'DisplayName', 'Monte Carlo Estimates');
hold on;

% Draw the latest CRB Ellipse
sigma_theta = sqrt(var_theta_crb);
sigma_phi   = sqrt(var_phi_crb);
t = linspace(0, 2*pi, 200);
ellipse_theta = rad2deg(sigma_factor * sigma_theta) * cos(t);
ellipse_phi   = rad2deg(sigma_factor * sigma_phi)   * sin(t);

plot(ellipse_theta, ellipse_phi, 'r-', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('CRB %d-\\sigma Region', sigma_factor));

xlabel('Theta Error (deg)');
ylabel('Phi Error (deg)');
title(sprintf('Error Distribution at SNR = %d dB', snr_db_sweep(end)));
legend('Location', 'best');
grid on; axis equal;
hold off;


%% LOCAL HELPER FUNCTIONS

function val = music_cost_function(ang, px, py, UU)
    % Evaluates the MUSIC pseudo-spectrum cost (to be minimized)
    theta = ang(1);
    phi   = ang(2);
    
    u = sin(theta) * cos(phi);
    v = sin(theta) * sin(phi);
    
    phase = 2 * pi * (px * u + py * v);
    a = exp(1j * phase);
    
    val = real(a' * UU * a);
end

function R_ss = spatial_smoothing(X, M, N, M_sub, N_sub)
    % Applies 2D Spatial Smoothing to the covariance matrix
    [~, N_snaps] = size(X);
    X_mat = reshape(X, [M, N, N_snaps]);
    
    Lx = M - M_sub + 1;
    Ly = N - N_sub + 1;
    L_total = Lx * Ly;
    
    R_ss = zeros(M_sub*N_sub, M_sub*N_sub);
    
    for i = 1:Lx
        for j = 1:Ly
            x_sub = X_mat(i:i+M_sub-1, j:j+N_sub-1, :);
            x_vec = reshape(x_sub, [M_sub*N_sub, N_snaps]);
            R_ss = R_ss + (x_vec * x_vec') / N_snaps;
        end
    end
    R_ss = R_ss / L_total;
end

function [var_theta, var_phi] = crb(theta_deg, phi_deg, snr_db, N_snapshots, px, py)
    % Computes the theoretical Cramer-Rao Bound for a single source
    M = length(px); 
    if length(py) ~= M
        error('Position vectors px and py must have the same length.');
    end
    
    theta = deg2rad(theta_deg);
    phi   = deg2rad(phi_deg);
    snr = 10^(snr_db/10);
    
    % Direction cosines
    u = sin(theta) * cos(phi);
    v = sin(theta) * sin(phi);
    
    % Steering vector
    phase = 2 * pi * (px * u + py * v);
    a = exp(1j * phase);
    
    % Derivatives
    du_dtheta = cos(theta) * cos(phi);
    dv_dtheta = cos(theta) * sin(phi);
    du_dphi = -sin(theta) * sin(phi);
    dv_dphi =  sin(theta) * cos(phi);
    
    % Derivative matrix D
    D = [1j*2*pi*(px*du_dtheta + py*dv_dtheta).*a, ...
         1j*2*pi*(px*du_dphi   + py*dv_dphi).*a];
         
    % Orthogonal projection matrix
    P = eye(M) - (a * a') / M;
    
    % Fisher Information Matrix (FIM) and CRB
    FIM = 2 * N_snapshots * snr * real(D' * P * D);
    CRB = pinv(FIM);
    
    var_theta = CRB(1,1);
    var_phi   = CRB(2,2);
end