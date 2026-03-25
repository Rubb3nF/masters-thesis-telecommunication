clear; clc; close all;

% This script simulates a target (e.g., an aircraft) flying in a straight 
% line and evaluates methods for predicting its future angular position 
% (Azimuth and Elevation) relative to an observer. 
% 
% It features two main sections:
%   1. A basic demonstration of the prediction using Rodrigues' rotation formula.
%   2. A quantitative comparison between a standard linear extrapolation 
%      and the Rodrigues-based spherical kinematics approach, including 
%      error analysis during the fastest part of the flyby.

%% Basic Scenario

% Scenario configuration 
v_plane = [200; 0; 0]; 
pos_plane_init = [-2000; 0; 1000];
dt = 0.5;
pos_target = [0; 200; 0]; 
N = 30;

true_theta = zeros(1, N);
true_phi   = zeros(1, N);
pred_theta = zeros(1, N);
pred_phi   = zeros(1, N);

pos_plane = pos_plane_init;

% Trajectory and predict generation 
for t = 1:N
    r_vec = pos_target - pos_plane;
    
    x = r_vec(1); 
    y = r_vec(2); 
    z = r_vec(3);
    dist = norm(r_vec);
    
    %  Radians
    true_theta(t) = acos(z/dist);     
    true_phi(t)   = atan2(y, x);      
    
    % PREDICTION 
    if t > 3
        % Assuming constant dt for dt1, dt2, dt_pred in this simple block
        [p_t, p_p] = predictNextPositionRodrigues( ...
            true_theta(t-3), true_phi(t-3), ...
            true_theta(t-2), true_phi(t-2), ...
            true_theta(t-1), true_phi(t-1), ...
            dt, dt, dt); 
        
        pred_theta(t) = p_t;
        pred_phi(t)   = p_p;
    else
        pred_theta(t) = nan; 
        pred_phi(t)   = nan;
    end
    
    pos_plane = pos_plane + v_plane * dt;
end

% Visualize Results 
time = (1:N) * dt;
figure('Name', 'Part 1: Rodrigues Prediction', 'Position', [100 100 1000 600]);

% Convert to degrees only for visualization
true_phi_deg   = rad2deg(true_phi);
pred_phi_deg   = rad2deg(pred_phi);
true_theta_deg = rad2deg(true_theta);
pred_theta_deg = rad2deg(pred_theta);

subplot(2,1,1);
plot(time, true_phi_deg, 'b.-', 'LineWidth', 1.5, 'DisplayName', 'True'); hold on;
plot(time, pred_phi_deg, 'rx--', 'LineWidth', 1, 'DisplayName', 'Prediction (Rodrigues)');
grid on;
title('Azimuth (Phi)');
ylabel('Degrees');
legend('Location', 'best');

subplot(2,1,2);
plot(time, true_theta_deg, 'b.-', 'LineWidth', 1.5, 'DisplayName', 'True'); hold on;
plot(time, pred_theta_deg, 'rx--', 'LineWidth', 1, 'DisplayName', 'Prediction (Rodrigues)');
grid on;
title('Elevation (Theta)');
ylabel('Degrees'); 
xlabel('Time (s)');

% Error in degrees (ignoring NaNs)
err_p = abs(true_phi_deg(4:end) - pred_phi_deg(4:end));
fprintf('Maximum error in Phi during the fast zone: %.4f degrees\n\n', max(err_p));


%% PART 2: Comparison (Linear vs. Rodrigues)

% Scenario Configuration
v_plane2 = [150; 0; 0];
pos_plane2 = [-1000; 100; 1000];
pos_target2 = [0; 0; 0];
dt2 = 1.0;
N2 = 25;

history_true = zeros(2, N2); 
history_lin  = zeros(2, N2); 
history_rod  = zeros(2, N2); 

% Simulation Loop
current_pos = pos_plane2;
for t = 1:N2
    
    % Ground Truth 
    r = pos_target2 - current_pos;
    dist = norm(r);
    
    th_true = acos(r(3)/dist);    
    phi_true = atan2(r(2), r(1));  
    
    history_true(:, t) = [th_true; phi_true];
    
    % Predictions
    if t > 3
        th_0 = history_true(1, t-3);
        ph_0 = history_true(2, t-3);
        
        th_1 = history_true(1, t-2);
        ph_1 = history_true(2, t-2);
        
        th_2 = history_true(1, t-1);
        ph_2 = history_true(2, t-1);
        
        % LINEAR PREDICTION
        d1_th = th_1 - th_0;
        d2_th = th_2 - th_1;
        
        d1_ph = ph_1 - ph_0;
        d2_ph = ph_2 - ph_1;
        
        beta = 0.7;
        
        pred_lin_th  = th_2 + d2_th + beta*(d2_th - d1_th);
        pred_lin_phi = ph_2 + d2_ph + beta*(d2_ph - d1_ph);
        
        % Wrap phi to [-pi, pi]
        pred_lin_phi = atan2(sin(pred_lin_phi), cos(pred_lin_phi));
        
        history_lin(:, t) = [pred_lin_th; pred_lin_phi];
        
        % RODRIGUES PREDICTION
        delta_t1 = dt2;      
        delta_t2 = dt2;      
        delta_t_pred = dt2;  
        
        [pred_rod_th, pred_rod_phi] = predictNextPositionRodrigues( ...
            th_0, ph_0, ...
            th_1, ph_1, ...
            th_2, ph_2, ...
            delta_t1, delta_t2, delta_t_pred);
        
        history_rod(:, t) = [pred_rod_th; pred_rod_phi];
        
    else
        history_lin(:, t) = [NaN; NaN];
        history_rod(:, t) = [NaN; NaN];
    end
    
    current_pos = current_pos + v_plane2 * dt2;
end

% Error Analysis
valid_idx = 4:N2;
true_phi_deg2 = rad2deg(history_true(2, valid_idx));
lin_phi_deg   = rad2deg(history_lin(2, valid_idx));
rod_phi_deg   = rad2deg(history_rod(2, valid_idx));

err_lin_phi = abs(true_phi_deg2 - lin_phi_deg);
err_rod_phi = abs(true_phi_deg2 - rod_phi_deg);

% Correction for 360-degree wrapping jumps
err_lin_phi(err_lin_phi > 300) = abs(err_lin_phi(err_lin_phi > 300) - 360);
err_rod_phi(err_rod_phi > 300) = abs(err_rod_phi(err_rod_phi > 300) - 360);

true_theta_deg2 = rad2deg(history_true(1, valid_idx));
lin_theta_deg   = rad2deg(history_lin(1, valid_idx));
rod_theta_deg   = rad2deg(history_rod(1, valid_idx));

err_lin_theta = abs(true_theta_deg2 - lin_theta_deg);
err_rod_theta = abs(true_theta_deg2 - rod_theta_deg);

fprintf('Mean Linear Error (Phi):    %.4f deg\n', mean(err_lin_phi));
fprintf('Mean Rodrigues Error (Phi): %.4f deg\n', mean(err_rod_phi));
fprintf('Mean Linear Error (Theta):    %.4f deg\n', mean(err_lin_theta));
fprintf('Mean Rodrigues Error (Theta): %.4f deg\n', mean(err_rod_theta));

% Visualization
figure('Name', 'Part 2: Linear vs Rodrigues', 'Position', [150 150 1200 800]);

% PHI
subplot(2, 2, 1);
plot(valid_idx, true_phi_deg2, 'k-', 'LineWidth', 2, 'DisplayName', 'True'); hold on;
plot(valid_idx, lin_phi_deg, 'b--o', 'DisplayName', 'Linear');
plot(valid_idx, rod_phi_deg, 'r-x', 'LineWidth', 1.5, 'DisplayName', 'Rodrigues');
grid on;
title('Azimuth Trajectory (Phi)');
ylabel('Degrees'); xlabel('Time step');
legend('Location', 'best');

subplot(2, 2, 2);
plot(valid_idx, err_lin_phi, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Linear Error'); hold on;
plot(valid_idx, err_rod_phi, 'r-x', 'LineWidth', 1.5, 'DisplayName', 'Rodrigues Error');
grid on;
title('Absolute Error in Phi');
ylabel('Error (Degrees)'); xlabel('Time step');
legend('Location', 'best');

% THETA
subplot(2, 2, 3);
plot(valid_idx, true_theta_deg2, 'k-', 'LineWidth', 2, 'DisplayName', 'True'); hold on;
plot(valid_idx, lin_theta_deg, 'b--o', 'DisplayName', 'Linear');
plot(valid_idx, rod_theta_deg, 'r-x', 'LineWidth', 1.5, 'DisplayName', 'Rodrigues');
grid on;
title('Elevation Trajectory (Theta)');
ylabel('Degrees'); xlabel('Time step');
legend('Location', 'best');

subplot(2, 2, 4);
plot(valid_idx, err_lin_theta, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Linear Error'); hold on;
plot(valid_idx, err_rod_theta, 'r-x', 'LineWidth', 1.5, 'DisplayName', 'Rodrigues Error');
grid on;
title('Absolute Error in Theta');
ylabel('Error (Degrees)'); xlabel('Time step');
legend('Location', 'best');


%% LOCAL FUNCTIONS

function [theta_next, phi_next] = predictNextPositionRodrigues( ...
    theta_prev2, phi_prev2, ...
    theta_prev1, phi_prev1, ...
    theta_curr,  phi_curr, ...
    dt1, dt2, dt_pred)

    % Description: Predicts the next angular position using Rodrigues' rotation
    %              formula based on the last three observations.-
    
    beta = 0.7;                 
    max_step = 45 * pi / 180;       
    eps_axis = 1e-8;
    
    % Conversion to unit vectors
    v0 = [sin(theta_prev2)*cos(phi_prev2);
          sin(theta_prev2)*sin(phi_prev2);
          cos(theta_prev2)];
          
    v1 = [sin(theta_prev1)*cos(phi_prev1);
          sin(theta_prev1)*sin(phi_prev1);
          cos(theta_prev1)];
          
    v2 = [sin(theta_curr)*cos(phi_curr);
          sin(theta_curr)*sin(phi_curr);
          cos(theta_curr)];
          
    % Rotation 0 -> 1
    k1 = cross(v0, v1);
    n1 = norm(k1);
    if n1 < eps_axis
        omega1 = [0; 0; 0];
    else
        k1 = k1 / n1;
        a1 = acos(max(min(dot(v0, v1), 1), -1));
        omega1 = a1 * k1;
    end
    
    % Rotation 1 -> 2
    k2 = cross(v1, v2);
    n2 = norm(k2);
    if n2 < eps_axis
        omega2 = [0; 0; 0];
    else
        k2 = k2 / n2;
        a2 = acos(max(min(dot(v1, v2), 1), -1));
        omega2 = a2 * k2;
    end
    
    % Convert to true angular velocity
    Omega1 = omega1 / dt1;
    Omega2 = omega2 / dt2;
    
    % Angular acceleration
    accel = (Omega2 - Omega1) / dt2;
    
    % Predict angular velocity
    Omega3 = Omega2 + beta * accel * dt_pred;
    
    % Revert to angular increment
    omega3 = Omega3 * dt_pred;
    alpha3 = norm(omega3);
    alpha3 = min(alpha3, max_step);
    
    if alpha3 < eps_axis
        v3 = v2;
    else
        k3 = omega3 / norm(omega3);
        % Rodrigues' rotation formula
        v3 = v2*cos(alpha3) + ...
             cross(k3, v2)*sin(alpha3) + ...
             k3*dot(k3, v2)*(1 - cos(alpha3));
        v3 = v3 / norm(v3);
    end
    
    % Return to spherical coordinates
    theta_next = acos(v3(3));
    phi_next   = atan2(v3(2), v3(1));
end