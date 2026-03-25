% This script simulates a 3D Direction of Arrival scenario. It generates 
% synthetic narrowband signals for multiple sources with varying SNRs.
% It then applies the 2D MUSIC algorithm to estimate the number of sources
% and their respective Zenith and Azimuth angles. Finally, it visualizes
% the 2D spatial pseudo-spectrum, marking the true and estimated DOAs.

clear; clc; close all;

% Parameters
M = 16;
N = 16;
fc = 1e9;
d = 0.5;
c = 3e8;

% Initialize Array
array = UniformPlanarArray(M, N, d, d, fc, c);
coherence = true;

% % Signal Model
signal_model = SignalModel(array);
K = 3; % Number of sources

% Ground Truth DOAs (Zenith, Azimuth)
thetas_true = deg2rad([70 30 40]); 
phis_true = deg2rad([-14 40 40]);

N_snapshots = 300; 
snr_db = [10 10 10]; % SNR per source
seed = 42;

fprintf('Generating signal for VCA with K=2 sources...\n');
[X, ~, ~] = signal_model.generateSignals(thetas_true, phis_true, N_snapshots, snr_db, ...
    'complex_sinusoid', 0.0, 'white', seed);

% Signal visualization
figure;
plot(fftshift(abs(fft(X(1,:)))));
title('Frequency Magnitude of the Received Signal (First Element)');
xlabel('Frequency Bin');
ylabel('Magnitude');

% Estimator Initialization
estimador = fast_MUSIC(array);
thetas_grid = linspace(0, pi/2, 91);
phis_grid = linspace(-pi, pi, 361);

% Source Enumeration
K_ER = estimador.calculateNumSources(X, coherence, true);
fprintf('%d sources detected \n', K_ER);

fprintf('Executing MUSIC estimation...\n');
tic;
[thetas_est, phis_est] = estimador.estimate(X, 3, coherence);

toc;

% Display Results
fprintf('True DOAs:\n');
for i = 1:K
    fprintf('  Source %d (Th, Ph): (%.2f°, %.2f°)\n', i, rad2deg(thetas_true(i)), rad2deg(phis_true(i)));
end

fprintf('Estimated DOAs:\n');
for i = 1:length(thetas_est)
    fprintf('  Peak %d (Th, Ph): (%.2f°, %.2f°)\n', i, rad2deg(thetas_est(i)), rad2deg(phis_est(i)));
end

% 2D Spectrum Calculation
spectrum_2d = estimador.musicSpectrumVectorized(X, 3, thetas_grid, phis_grid, coherence);
spectrum_dB = 10 * log10(spectrum_2d / max(spectrum_2d(:)));

thetas_grid_deg = rad2deg(thetas_grid);
phis_grid_deg = rad2deg(phis_grid);
[THETA_GRID_DEG, PHI_GRID_DEG] = meshgrid(thetas_grid_deg, phis_grid_deg);

% Visualization
figure;
surf(THETA_GRID_DEG, PHI_GRID_DEG, spectrum_dB, 'EdgeColor', 'none');
colormap jet; 
title('2D MUSIC Spectrum (Zenith vs. Azimuth)');
xlabel('Zenith Angle \theta (degrees)');
ylabel('Azimuth Angle \phi (degrees)');
zlabel('MUSIC Power (dB)');
colorbar;
view(45, 30); 
axis tight;
hold on;

% True DOAs
for k = 1:K
    [~, theta_idx] = min(abs(thetas_grid - thetas_true(k)));
    [~, phi_idx]   = min(abs(phis_grid - phis_true(k)));
    z_true = spectrum_dB(phi_idx, theta_idx);
    plot3(rad2deg(thetas_true(k)), rad2deg(phis_true(k)), z_true, ...
          'wo', 'MarkerSize', 11, 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'True DOA');
end

% Estimated DOAs
K_est = length(thetas_est);
for k = 1:K_est
    theta_est = thetas_est(k);
    phi_est   = phis_est(k);  
    [~, theta_idx_est] = min(abs(thetas_grid - theta_est));
    [~, phi_idx_est]   = min(abs(phis_grid - phi_est));
    z_est = spectrum_dB(phi_idx_est, theta_idx_est);
    plot3(rad2deg(theta_est), rad2deg(phi_est), z_est, ...
          'ko', 'MarkerSize', 11, 'LineWidth', 2, 'MarkerFaceColor', 'g', 'DisplayName', 'Estimated DOA');
end
hold off;