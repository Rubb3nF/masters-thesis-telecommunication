% Evaluates the Probability of Detection (Pd) for the number of sources
% using Eigenvalue Ratio (ER), Akaike Information Criterion (AIC), and 
% Minimum Description Length (MDL) across a sweep of SNR values.

clear; clc; close all;
%% System & Array Configuration 

M = 16;
N = 16; % Total antennas = 256
fc = 1e9;
d = 0.5;
c = 3e8;

% Initialize Array and Signal Model
array_vca = UniformPlanarArray(M, N, d, d, fc, c);
signal_model = SignalModel(array_vca);

% Simulation Parameters
K_true = 2; % True number of sources
thetas_true = deg2rad([70 30]); 
phis_true   = deg2rad([-14 40]);
N_snapshots = 300; 

%% Monte Carlo Setup
snr_range = -25:2:10; % SNR sweep from -25dB to 10dB
n_montecarlo = 100;   % Iterations per SNR point

% Arrays to store probability of detection (0 to 1)
prob_ER  = zeros(length(snr_range), 1);
prob_AIC = zeros(length(snr_range), 1);
prob_MDL = zeros(length(snr_range), 1);


%% Main Loop 
for i = 1:length(snr_range)
    current_snr = snr_range(i);
    snr_vec = [current_snr, current_snr + 10]; % Same base SNR for both sources
    
    % Hit counters
    hits_ER = 0;
    hits_AIC = 0;
    hits_MDL = 0;
    
    for mc = 1:n_montecarlo
        % Generate random seed so noise changes in each iteration
        seed = randi(100000); 
        
        % Generate signal 
        [X, ~, ~] = signal_model.generateSignals(thetas_true, phis_true, N_snapshots, snr_vec, ...
            'complex_sinusoid', 0.0, 'white', seed);
        
        % EIGENVALUE COMPUTATION
        [n_sensors, n_snaps] = size(X);
        Rxx = (X * X') / n_snaps;
        [~, D] = eig(Rxx);
        lambdas = sort(diag(D), 'descend'); % Sort from largest to smallest
        
        % METHOD 1: Eigenvalue Ratio (ER) 
        % Avoid division by zero or NaN in extreme cases
        lambdas_safe = lambdas; 
        lambdas_safe(lambdas_safe < 1e-9) = 1e-9; 
        ratios = lambdas_safe(1:end-1) ./ lambdas_safe(2:end);
        [~, k_er_idx] = max(ratios);
        K_est_ER = k_er_idx; % The index of the maximum jump is K
        
        % PREPARATION FOR AIC and MDL 
        aic_vals = zeros(n_sensors-1, 1);
        mdl_vals = zeros(n_sensors-1, 1);
        
        for k = 0:n_sensors-2
            % M is the total number of sensors, m is (M-k)
            m = n_sensors - k; 
            
            % Eigenvalues of the noise subspace (lambda_{k+1} ... lambda_M)
            noise_eigs = lambdas(k+1:end);
            
            % Numerator: Geometric Mean (Product^(1/m))
            numerador = prod(noise_eigs)^(1/m);
            
            % Denominator: Arithmetic Mean (Sum / m)
            denominador = mean(noise_eigs);
            
            % Ratio (always <= 1 due to AM-GM inequality)
            ratio = numerador / denominador;
            
            % Protection against log(0) if the ratio is extremely small
            if ratio <= 0
                log_ratio = -1e10; 
            else
                log_ratio = log(ratio);
            end
            
            % Calculate the Log-Likelihood Term
            % ((M-k)*N) * log(ratio)
            log_likelihood_term = (m * n_snaps) * log_ratio;
            
            % AIC Criterion
            % Formula: -2 * log(...) + 2k(2M-k)
            penalty_aic = 2 * k * (2*n_sensors - k);
            aic_vals(k+1) = -2 * log_likelihood_term + penalty_aic;
            
            % MDL Criterion 
            % -1 * log(...) + 0.5 * k * (2M-k) * log(N)
            penalty_mdl = 0.5 * k * (2*n_sensors - k) * log(n_snaps);
            mdl_vals(k+1) = -1 * log_likelihood_term + penalty_mdl;
        end
        
        [~, idx_aic] = min(aic_vals);
        K_est_AIC = idx_aic - 1; % Subtract 1 because index 1 corresponds to k=0
        
        [~, idx_mdl] = min(mdl_vals);
        K_est_MDL = idx_mdl - 1;
        
        % HIT VERIFICATION 
        if K_est_ER == K_true,  hits_ER = hits_ER + 1;  end
        if K_est_AIC == K_true, hits_AIC = hits_AIC + 1; end
        if K_est_MDL == K_true, hits_MDL = hits_MDL + 1; end
    end
    
    % Calculate percentages
    prob_ER(i)  = hits_ER / n_montecarlo;
    prob_AIC(i) = hits_AIC / n_montecarlo;
    prob_MDL(i) = hits_MDL / n_montecarlo;
    
    fprintf('SNR: %3d dB completed. Hits -> ER: %.2f, AIC: %.2f, MDL: %.2f\n', ...
        current_snr, prob_ER(i), prob_AIC(i), prob_MDL(i));
end

%% Visualization
figure('Color', 'w', 'Name', 'Source Estimation Performance');
plot(snr_range, prob_ER, '-o', 'LineWidth', 2, 'DisplayName', 'Eigenvalue Ratio (ER)');
hold on;
plot(snr_range, prob_AIC, '-^', 'LineWidth', 2, 'DisplayName', 'AIC');
plot(snr_range, prob_MDL, '-s', 'LineWidth', 2, 'DisplayName', 'MDL');
hold off;

grid on;
xlabel('SNR (dB)');
ylabel('Detection Probability (Pd)');
title(sprintf('Performance Comparison of Source Estimation Methods (N_{snap}=%d, K=%d)', N_snapshots, K_true));
legend('Location', 'best');
ylim([-0.05 1.05]);