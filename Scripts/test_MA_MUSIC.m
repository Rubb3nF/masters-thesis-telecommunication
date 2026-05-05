% Script: Dynamic DOA Tracking Scenarios (MUSIC)
%
% Description:
% Evaluates the tracking performance of the motion_MUSIC estimator under 
% three different dynamic scenarios:
%   1. 'APPEARANCE' : A 3rd source suddenly appears at step 15.
%   2. 'OVERLAP'    : Sources start from the same origin and diverge.
%   3. 'IDEAL'      : Non-intersecting, well-separated trajectories.

clear; clc; close all;
%% Scenario Selection & General Configuration
% Choose scenario: 1 ('APPEARANCE'), 2 ('OVERLAP'), or 3 ('IDEAL')
scenario_mode = 3; 

% Array & System Setup
M = 16; N = 16; fc = 1e9; d = 0.5; c = 3e8;

% Dynamic Parameters
K_max = 3;              
N_snapshots = 500;
snr_db = 20;        
coherent = true;

% Different number of steps depending on scenario
if scenario_mode == 1
    num_steps = 40;
else
    num_steps = 45;
end

%% Object Instantiation
array_vca = UniformPlanarArray(M, N, d, d, fc, c);
signal_model = SignalModel(array_vca);

estimator = motion_MUSIC_opt(array_vca);
estimator.history = {};   

%% Trajectory Generation (Scenario Switch) 
path_theta_deg = zeros(num_steps, K_max);
path_phi_deg   = zeros(num_steps, K_max);
active_sources = 3 * ones(num_steps, 1); % Tracks how many sources exist at step k

switch scenario_mode
    case 1 % SCENARIO 1: DYNAMIC APPEARANCE 
        aparece_en_paso = 15;
        active_sources(1:aparece_en_paso-1) = 2; % Only 2 sources initially
        t_vec = linspace(0, 2*pi, num_steps);
        
        path_theta_deg(:,1) = 30 + 10*sin(t_vec);
        path_phi_deg(:,1)   = -40 + 30*cos(t_vec);
        path_theta_deg(:,2) = linspace(60, 20, num_steps);
        path_phi_deg(:,2)   = linspace(20, -20, num_steps);
        path_theta_deg(:,3) = linspace(10, 50, num_steps); % Surprise source
        path_phi_deg(:,3)   = linspace(80, 40, num_steps);

    case 2 % SCENARIO 2: STARTING OVERLAP
        t_vec = linspace(0, 1, num_steps);
        theta_start = 30; phi_start = 0;
        
        % Source 1: Top-Left
        path_theta_deg(:,1) = linspace(theta_start, 15, num_steps); 
        path_phi_deg(:,1)   = linspace(phi_start, -50, num_steps);
        % Source 2: Independent (Bottom-Right)
        path_theta_deg(:,2) = 65 + 10*cos(linspace(0, 2*pi, num_steps));
        path_phi_deg(:,2)   = 50 + 15*sin(linspace(0, 2*pi, num_steps));
        % Source 3: Vertical movement from origin
        path_theta_deg(:,3) = linspace(theta_start, 75, num_steps); 
        path_phi_deg(:,3)   = linspace(phi_start, 10, num_steps); 

    case 3 % SCENARIO 3: IDEAL TRAJECTORIES
        t_vec = linspace(0, 2*pi, num_steps);
        
        % Source 1: Top-Left
        path_theta_deg(:,1) = 25 + 8*sin(t_vec);
        path_phi_deg(:,1)   = -60 + 10*cos(t_vec);
        % Source 2: Bottom-Right
        path_theta_deg(:,2) = 65 + 5*cos(t_vec);
        path_phi_deg(:,2)   = 50 + 15*sin(t_vec);
        % Source 3: Central
        path_theta_deg(:,3) = linspace(30, 60, num_steps); 
        path_phi_deg(:,3)   = linspace(0, 5, num_steps);
end

path_theta = deg2rad(path_theta_deg);
path_phi   = deg2rad(path_phi_deg);

%% Initialization of Logging Metrics 
log_theta_est   = nan(num_steps, K_max);
log_phi_est     = nan(num_steps, K_max);
log_total_error = nan(num_steps, K_max);
execution_times = zeros(num_steps, 1);

% Visualization Grids
grid_step_vis = 1.0; 
thetas_vis_deg = 0:grid_step_vis:90;
phis_vis_deg   = -180:grid_step_vis:180;
thetas_vis = deg2rad(thetas_vis_deg);
phis_vis   = deg2rad(phis_vis_deg);

figure('Color', 'white', 'Position', [100, 100, 1100, 600], 'Name', 'Real-Time Tracking');

%% Main Tracking Loop
for k = 1:num_steps
    K_actual = active_sources(k);
    
    current_th = path_theta(k, 1:K_actual);
    current_ph = path_phi(k, 1:K_actual);
    
    % Generate signals
    [X, ~, ~] = signal_model.generateSignals(current_th, current_ph, ...
        N_snapshots, snr_db, 'complex_sinusoid', 0, 'white', k);
        
    % LATENCY MEASUREMENT (Estimation Only)
    % Source number estimation
    if ismethod(estimator, 'calculateNumSources')
        K_est = estimator.calculateNumSources(X, true, false); 
    else
        K_est = K_actual;
    end
    
    tic;
    [th_est, ph_est] = estimator.estimate(X, K_est, 0.1*pi/180, coherent);
    execution_times(k) = toc; 
    
    % LOGGING & DATA ASSOCIATION (Euclidean Distance) 
    if ~isempty(th_est)
        estimates = [th_est(:), ph_est(:)];
        truths    = [current_th(:), current_ph(:)];
        dist_mat = pdist2(estimates, truths);
        
        for i = 1:K_actual
            if all(isinf(dist_mat(:))), break; end
            [min_val, min_idx] = min(dist_mat(:));
            [r, c] = ind2sub(size(dist_mat), min_idx);
            
            log_theta_est(k, c) = rad2deg(th_est(r));
            log_phi_est(k, c)   = rad2deg(ph_est(r));
            log_total_error(k, c) = rad2deg(min_val); % Total Euclidean error in degrees
            
            % Blackout used rows/columns
            dist_mat(r, :) = Inf; dist_mat(:, c) = Inf;
        end
    end
    
    % REAL-TIME VISUALIZATION
    clf; hold on;
    try
        spec_full = estimator.musicSpectrumVectorized(X, K_actual, thetas_vis, phis_vis, coherent);
        spec_db = 10*log10(spec_full / max(spec_full(:)));
        imagesc(phis_vis_deg, thetas_vis_deg, spec_db'); 
        axis xy; colormap jet; caxis([-20 0]); colorbar;
    catch
        axis([-100 100 0 90]);
    end
    
    % Draw true positions
    plot(path_phi_deg(k, 1), path_theta_deg(k, 1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'True 1');
    plot(path_phi_deg(k, 2), path_theta_deg(k, 2), 'bs', 'MarkerFaceColor', 'b', 'DisplayName', 'True 2');
    if K_actual == 3
        plot(path_phi_deg(k, 3), path_theta_deg(k, 3), 'ms', 'MarkerFaceColor', 'm', 'DisplayName', 'True 3');
    end
    
    % Draw estimates
    if ~isempty(th_est)
        plot(rad2deg(ph_est), rad2deg(th_est), 'wx', 'MarkerSize', 12, 'LineWidth', 2.5, 'DisplayName', 'Estimations');
    end
    
    % Draw tracker ROI ellipses (if applicable)
    if isprop(estimator, 'last_pred_theta') && ~isempty(estimator.last_pred_theta)
        ang = linspace(0, 2*pi, 50);
        for src = 1:length(estimator.last_pred_theta)
            c_th = rad2deg(estimator.last_pred_theta(src));
            c_ph = rad2deg(estimator.last_pred_phi(src));
            r_th = rad2deg(estimator.last_r_theta(src));
            r_ph = rad2deg(estimator.last_r_phi(src));
            h_roi = plot(c_ph + r_ph*cos(ang), c_th + r_th*sin(ang), 'r-', 'LineWidth', 2);
            if src == 1, set(h_roi, 'DisplayName', 'ROI'); else, set(get(get(h_roi,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); end
        end
    end
    
    legend('Location', 'northeastoutside', 'TextColor', 'black', 'FontSize', 10);
    title(sprintf('Step %d/%d | Latency: %.2f ms', k, num_steps, execution_times(k)*1000));
    xlabel('\phi Azimuth (deg)'); ylabel('\theta Zenith (deg)');
    xlim([-100 100]); ylim([0 90]); grid on;
    drawnow;
end

%% STATIC REPORT: TRAJECTORY & PERFORMANCE 

% Trajectory Reconstruction
figure('Color', 'white', 'Name', 'Trajectory Reconstruction', 'Position', [100 100 800 500]);
hold on;
colors = ['g', 'b', 'm'];
for src = 1:K_max
    plot(path_phi_deg(:, src), path_theta_deg(:, src), '--', 'Color', colors(src), 'LineWidth', 1.5, 'DisplayName', sprintf('True %d', src));
    plot(log_phi_est(:, src), log_theta_est(:, src), 'o', 'MarkerEdgeColor', colors(src), 'MarkerSize', 5, 'DisplayName', sprintf('Est %d', src));
end
xlabel('\phi Azimuth (deg)'); ylabel('\theta Zenith (deg)');
title('Reconstructed vs True Trajectories');
grid on; legend('Location', 'northeastoutside');

% Error & Latency Metrics
figure('Color', 'white', 'Name', 'Performance Metrics', 'Position', [150 150 800 600]);

subplot(2,1,1);
plot(log_total_error, 'LineWidth', 1.5);
yline(nanmean(log_total_error(:)), 'r--', 'Overall Mean Error', 'LineWidth', 1.5);
title('Total Angular Error (Euclidean Distance)');
ylabel('Error (deg)'); xlabel('Time Step');
grid on; legend('Source 1', 'Source 2', 'Source 3', 'Location', 'northeast');

subplot(2,1,2);
bar(execution_times * 1000, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'none');
hold on;
yline(mean(execution_times)*1000, 'r--', sprintf('Avg: %.2f ms', mean(execution_times)*1000), 'LineWidth', 2);
title('Algorithm Computational Latency (.estimate only)');
ylabel('Time (ms)'); xlabel('Time Step');
grid on;

% Final Console Print
fprintf('\n FINAL PERFORMANCE REPORT \n');
for src = 1:K_max
    rmse_th = sqrt(nanmean((log_theta_est(:, src) - path_theta_deg(:, src)).^2));
    fprintf('Source %d -> RMSE Theta: %.4f deg\n', src, rmse_th);
end
fprintf('Mean RMSE Total (All Sources): %.4f deg\n', nanmean(log_total_error(:)));
fprintf('Average Latency: %.2f ms per frame\n', mean(execution_times)*1000);
fprintf('Max Latency:     %.2f ms\n', max(execution_times)*1000);