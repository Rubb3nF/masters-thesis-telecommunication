% This script validates the Cramer-Rao Bound for Direction of Arrival 
% estimation using a Uniform Planar Array. It runs a Monte Carlo 
% simulation to estimate the DOA (Theta and Phi) of a single target using 
% 2D Spatial Smoothing and MUSIC. The resulting estimation errors are then 
% compared against the theoretical CRB confidence ellipse.

clear; clc; close all;

%% System Parameters

% Array and Signal Parameters
M = 16;             % Number of elements in x-direction
N = 16;             % Number of elements in y-direction
d = 0.5;            % Element spacing (wavelengths)
fc = 1e9;           % Carrier frequency (Hz)
c = 3e8;            % Speed of light (m/s)

% Spatial Smoothing Sub-array Dimensions
M_sub = 16;         
N_sub = 16;         

% Initialize Array and Signal Models
array_vca = UniformPlanarArray(M, N, d, d, fc, c);
signal_model = SignalModel(array_vca);
estimador = motion_MUSIC_opt(array_vca);

% Retrieve Sensor Positions
posiciones = array_vca.element_positions;
x_pos = posiciones(1, :)';
y_pos = posiciones(2, :)';

% Target True Position
theta_true_deg = 69.99;
phi_true_deg   = -14.04;
theta_true = deg2rad(theta_true_deg);
phi_true   = deg2rad(phi_true_deg);

% Simulation Parameters
N_snapshots = 300;
snr_db = -15;
num_montecarlo = 500;

% Storage for Errors
errors_theta = zeros(num_montecarlo, 1);
errors_phi   = zeros(num_montecarlo, 1);


%% Sub-array Geometry Configuration
vec_x_sub = (0:M_sub-1) - (M_sub-1)/2;
vec_y_sub = (0:N_sub-1) - (N_sub-1)/2;
[X_sub, Y_sub] = ndgrid(vec_x_sub, vec_y_sub);
p_x_sub = X_sub(:) * d;
p_y_sub = Y_sub(:) * d;


%% Monte Carlo Simulation Loop
wb = waitbar(0, 'Running Monte Carlo simulation');

for mc = 1:num_montecarlo
    waitbar(mc/num_montecarlo, wb);
    
    seed = 100 + mc; % Unique seed per iteration
    
    % Generate signal snapshots
    [X, ~, ~] = signal_model.generateSignals( ...
        theta_true, phi_true, ...
        N_snapshots, snr_db, ...
        'complex_sinusoid', 0.0, 'white', seed);
    
    % Initial Guess (MUSIC on full array)
    theta_grid = linspace(theta_true - 0.1, theta_true + 0.1, 20);
    phi_grid   = linspace(phi_true - 0.1, phi_true + 0.1, 20);
    
    try
        % Only one source is assumed
        K_ER = 1; 
        [theta_init, phi_init] = estimador.estimate(X, K_ER, theta_grid, phi_grid, true);
        start_point = [theta_init, phi_init];
    catch
        % Fallback if grid search fails
        start_point = [theta_true, phi_true];
    end
    
    % Refinement
    R_ss = spatial_smoothing(X, M, N, M_sub, N_sub);
    
    % Eigen-decomposition
    [E, Dvals] = eig(R_ss);
    [~, idx] = sort(diag(Dvals), 'ascend');
    
    % Noise subspace
    En = E(:, idx(1:end-1)); 
    UU = En * En';
    
    % Cost function for fminsearch
    cost_func = @(ang) music_cost_function(ang, p_x_sub, p_y_sub, UU);
    options = optimset('Display', 'off', 'TolX', 1e-9);
    
    % Refine DOA estimation
    refined_pos = fminsearch(cost_func, start_point, options);
    
    % Store absolute errors in radians
    errors_theta(mc) = refined_pos(1) - theta_true;
    errors_phi(mc)   = refined_pos(2) - phi_true;
end
close(wb);


%% CRB Computation
[var_theta_crb, var_phi_crb] = crb( ...
    theta_true_deg, phi_true_deg, snr_db, ...
    N_snapshots, x_pos, y_pos);

sigma_theta = sqrt(var_theta_crb);
sigma_phi   = sqrt(var_phi_crb);


%% Confidence Region Test
sigma_factor = 1; % 1-sigma ellipse
inside = (errors_theta.^2 ./ (sigma_factor^2 * var_theta_crb) + ...
          errors_phi.^2   ./ (sigma_factor^2 * var_phi_crb)) <= 1;
          
percentage_inside = 100 * sum(inside) / num_montecarlo;

fprintf('\n--- Simulation Results ---\n');
fprintf('True direction: Theta = %.2f deg, Phi = %.2f deg\n', theta_true_deg, phi_true_deg);
fprintf('Monte Carlo points inside 1-sigma CRB region: %.2f%%\n', percentage_inside);


%% Visualization
figure('Color', 'w', 'Name', 'CRB Confidence Validation', 'Position', [100 100 600 500]);
plot(rad2deg(errors_theta), rad2deg(errors_phi), 'bx', 'LineWidth', 1.2, 'DisplayName', 'Monte Carlo Estimates');
hold on;

% CRB Ellipse
t = linspace(0, 2*pi, 200);
ellipse_theta = rad2deg(1*sigma_theta) * cos(t);
ellipse_phi   = rad2deg(1*sigma_phi)   * sin(t);
plot(ellipse_theta, ellipse_phi, 'r-', 'LineWidth', 2.5, 'DisplayName', '1-\sigma CRB Region');

grid on; axis equal;
xlabel('Theta Error (deg)');
ylabel('Phi Error (deg)');
title(sprintf('CRB Confidence Region (%.2f%% Inside)', percentage_inside));
legend('Location', 'best');
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
            % Extract sub-array data
            x_sub = X_mat(i:i+M_sub-1, j:j+N_sub-1, :);
            x_vec = reshape(x_sub, [M_sub*N_sub, N_snaps]);
            
            % Accumulate covariance
            R_ss = R_ss + (x_vec * x_vec') / N_snaps;
        end
    end
    % Average over all sub-arrays
    R_ss = R_ss / L_total;
end

function [var_theta, var_phi] = crb(theta_deg, phi_deg, snr_db, N_snapshots, px, py)
    % Computes the theoretical Cramer-Rao Bound for a single source
    % px, py: Column vectors with sensor coordinates (in wavelengths)
    
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
    
    % Fisher Information Matrix and CRB
    FIM = 2 * N_snapshots * snr * real(D' * P * D);
    CRB = pinv(FIM);
    
    var_theta = CRB(1,1);
    var_phi   = CRB(2,2);
end