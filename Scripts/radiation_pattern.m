% This script compares a standard Uniform Planar Array (UPA) with a 
% Jittered UPA. It visualizes the physical placement of the antenna 
% elements and computes the 1D Array Factor (E-Plane cut) to demonstrate 
% how spatial jittering helps suppress grating lobes.
%
% The 1D Array Factor is calculated along the u-axis as:
% $$AF(u) = \sum_{n=1}^{N_{elements}} e^{j k x_n u}$$
% where $u = \sin(\theta)$ and $k = 2\pi/\lambda$.

clear; clc; close all;

% Parameters
M = 8; 
N = 8; 
dx = 1; 
dy = 1; 
jitter_val = 0.35; 
fc = 1e9; 
c = 3e8; 
wavelength = c / fc;
k = 2 * pi / wavelength;

% Generate Arrays 
% Standard UPA (jitter = 0)
upa_std = JitteredUPA(M, N, dx, dy, 0, fc, c);
% Jittered UPA
upa_jit = JitteredUPA(M, N, dx, dy, jitter_val, fc, c);

% Physical Geometry
figure('Color', 'w', 'Name', 'Antenna Geometry', 'Position', [100 100 500 400]);
plot(upa_std.element_positions(1,:), upa_std.element_positions(2,:), ...
    'rs', 'MarkerSize', 8, 'DisplayName', 'Original Grid'); 
hold on;
plot(upa_jit.element_positions(1,:), upa_jit.element_positions(2,:), ...
    'b.', 'MarkerSize', 12, 'DisplayName', 'Jittered Positions');
axis equal; 
grid on; 
xlabel('x (m)'); 
ylabel('y (m)');
title('Physical Antenna Placement'); 
legend('Location', 'best');
hold off;

% 1D Cut in Electrical Space (u-axis)
% u-axis definition (u = sin(theta))
u = linspace(-1.2, 1.2, 1000);

% Array Factor calculation (AF)
AF_std = abs(sum(exp(1j * k * (upa_std.element_positions(1,:)' * u)), 1));
AF_jit = abs(sum(exp(1j * k * (upa_jit.element_positions(1,:)' * u)), 1));

% dB normalization
AF_std_dB = 20 * log10(AF_std / max(AF_std));
AF_jit_dB = 20 * log10(AF_jit / max(AF_jit));

figure('Color', 'w', 'Name', 'E-Plane Cut', 'Position', [600 100 700 400]);
plot(u, AF_std_dB, 'r', 'LineWidth', 1.5, 'DisplayName', 'Standard UPA (Grating Lobes at 0dB)'); 
hold on;
plot(u, AF_jit_dB, 'b', 'LineWidth', 1.5, 'DisplayName', 'Jittered UPA (Suppressed)');

xline([-1, 1], '--k', 'Visible Region', 'LabelVerticalAlignment', 'bottom');

grid on; 
ylim([-30 0]); 
xlabel('u (\sin \theta)'); 
ylabel('Magnitude (dB)');
title('E-Plane Cut: Grating Lobe Suppression');
legend('Location', 'southoutside');
hold off;