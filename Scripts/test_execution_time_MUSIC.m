% Evaluates the computational cost (execution time) of the vectorized 
% MUSIC spectrum estimation across different Uniform Planar Array
% sizes and angular search grid resolutions.

clear; clc; close all;
%% General Parameters
fc = 1e9;
d  = 0.5;
c  = 3e8;

% Source Configuration
K_true = 2;
thetas_true = deg2rad([70 30]); 
phis_true   = deg2rad([-13 40]);
N_snapshots = 300;
snr_db = [10 10];
seed = 42;
coherence = false;


%% Benchmark Configuration 
% Array sizes to evaluate (M x N)
array_sizes = [4 8 16 32];

% Grid resolutions in degrees
grid_steps_deg = [5 3 2 1 0.75 0.5 0.3];

% Pre-allocate results matrix
execution_time = zeros(length(array_sizes), length(grid_steps_deg));


%% Main Benchmark Loop
for a = 1:length(array_sizes)
    M = array_sizes(a);
    N = array_sizes(a);
    
    % Initialize Array and Signal Model
    array_vca = UniformPlanarArray(M, N, d, d, fc, c);
    signal_model = SignalModel(array_vca);
    
    % Generate Signals
    [X, ~, ~] = signal_model.generateSignals( ...
        thetas_true, phis_true, N_snapshots, snr_db, ...
        'complex_sinusoid', 0.0, 'white', seed);
        
    % Initialize MUSIC Estimator
    estimator = MUSIC(array_vca);
    
    for g = 1:length(grid_steps_deg)
        step_deg = grid_steps_deg(g);
        
        % Angular Grids
        thetas_grid = deg2rad(0:step_deg:90);
        phis_grid   = deg2rad(-180:step_deg:180);
        
        % Execution Time Measurement
        tic;
        estimator.musicSpectrumVectorized(X, K_true, thetas_grid, phis_grid, coherence);
        execution_time(a, g) = toc;
    end
end

%% Plot Results
figure('Color', 'w', 'Name', 'MUSIC Execution Time Benchmark', 'Position', [150 150 700 500]);
hold on; grid on;

markers = {'o', 's', '^', 'd'};
colors = lines(length(array_sizes)); % Use distinct colors

for a = 1:length(array_sizes)
    plot(grid_steps_deg, execution_time(a,:), ...
        'LineWidth', 2, ...
        'Color', colors(a,:), ...
        'Marker', markers{a}, ...
        'MarkerFaceColor', colors(a,:), ...
        'DisplayName', sprintf('Array %dx%d', array_sizes(a), array_sizes(a)));
end

set(gca, 'XDir', 'reverse');
set(gca, 'YScale', 'log');  
xlabel('Search Grid Resolution (degrees)');
ylabel('Execution Time (s) [Log Scale]');
title('MUSIC Execution Time vs. Grid Resolution & Array Size');
legend('Location', 'northwest');
hold off;