% This script verifies the 'generateTrajectory' method of the 'pseudo_KF' 
% class. It generates and visualizes three distinct 3D flight paths 
% (linear, circular, and curved via splines) around a stationary target.
% It demonstrates how to parameterize and plot different observer/target 
% movements for tracking simulations.

clear; clc; close all;

% Create an instance of the class
kf_gen = pseudo_KF();

% Common Definitions 
N_SAMPLES = 200;
Target_Position = [0, 0, 0];

% Linear Trajectory
start_linear = [50, 10, 5];
end_linear   = [100, 20, 15];

traj_linear = kf_gen.generateTrajectory(start_linear, end_linear, ...
    'Mode', 'linear', ...
    'NSamples', N_SAMPLES);

% Circular Trajectory
start_circle = [50, 50, 30];
end_circle   = [50, 50, 30]; % The end position is ignored in 'circle' mode
circle_params = struct('radius', 40, 'angle', 270, 'direction', 'right');

traj_circle = kf_gen.generateTrajectory(start_circle, end_circle, ...
    'Mode', 'circle', ...
    'NSamples', N_SAMPLES, ...
    'Params', circle_params);

% Curved Trajectory (Spline)
start_curve = [-50, -50, 40];
end_curve   = [-10, 10, 10];
waypoints = [
    -60,   0, 50;  
    -30,  50, 25; 
      0,  20,  5     
];

traj_curve = kf_gen.generateTrajectory(start_curve, end_curve, ...
    'Mode', 'curve', ...
    'NSamples', N_SAMPLES, ...
    'waypoints', waypoints);


figure('Color', 'w', 'Name', 'Trajectory Generation Test');
hold on;
grid on;

% Target
plot3(Target_Position(1), Target_Position(2), Target_Position(3), ...
    'r*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Target (0,0,0)');

% Linear Trajectory
plot3(traj_linear(:, 1), traj_linear(:, 2), traj_linear(:, 3), ...
    'b-', 'LineWidth', 2, 'DisplayName', 'Linear');
plot3(start_linear(1), start_linear(2), start_linear(3), ...
    'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot3(end_linear(1), end_linear(2), end_linear(3), ...
    'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'HandleVisibility', 'off');

% Circular Trajectory
color_circle = [0.85, 0.325, 0.098]; % Standard MATLAB orange
plot3(traj_circle(:, 1), traj_circle(:, 2), traj_circle(:, 3), ...
    'Color', color_circle, 'LineWidth', 2, 'DisplayName', 'Circle');
plot3(start_circle(1), start_circle(2), start_circle(3), ...
    'o', 'MarkerFaceColor', color_circle, 'MarkerEdgeColor', color_circle, ...
    'MarkerSize', 6, 'HandleVisibility', 'off');

% Curved Trajectory
plot3(traj_curve(:, 1), traj_curve(:, 2), traj_curve(:, 3), ...
    'g-', 'LineWidth', 2, 'DisplayName', 'Curve (Spline)');
plot3(start_curve(1), start_curve(2), start_curve(3), ...
    'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot3(end_curve(1), end_curve(2), end_curve(3), ...
    'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'HandleVisibility', 'off');
plot3(waypoints(:,1), waypoints(:,2), waypoints(:,3), ...
    'kx', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Waypoints');

title('Trajectory Comparison (Linear, Circle, Curve)');
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
zlabel('Z-axis (m)');
legend('Location', 'best');
view(3);
hold off;