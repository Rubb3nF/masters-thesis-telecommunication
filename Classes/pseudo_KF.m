classdef pseudo_KF
    % Description: Target tracking class implementing multiple Kalman Filter 
    %              variants (PLKF, WUPLKF, UKF) and trajectory generation.

    
    properties
    end
    
    methods
        function obj = pseudo_KF()
            % Constructor
        end
        
        function trajectory = generateTrajectory(obj, startPos, endPos, varargin)

            % Description: Generates a 3D trajectory for the target or observer.
            %
            % Inputs:
            %   startPos - Initial 3D position [x, y, z]
            %   endPos   - Final 3D position [x, y, z]
            %   varargin - Key-value pairs for trajectory parameters:
            %              'Mode': 'linear', 'circle', or 'curve'
            %              'waypoints': Matrix of intermediate points
            %              'Duration': Total duration (seconds)
            %              'NSamples': Number of generated points
            %              'Params': Struct with maneuver parameters
            %              'Speed', 'SpeedProfile', 'noise_level'
            %
            % Outputs:
            %   trajectory - Generated Nx3 matrix of [x, y, z] coordinates

            p = inputParser;
            addParameter(p, 'Mode', 'linear'); 
            addParameter(p, 'waypoints', []);
            addParameter(p, 'Duration', 60); 
            addParameter(p, 'NSamples', 1000); 
            addParameter(p, 'Params', struct()); 
            addParameter(p, 'Speed', 200);
            addParameter(p, 'SpeedProfile', []);
            addParameter(p, 'noise_level', 0);
            parse(p, varargin{:});
            
            switch p.Results.Mode
                case 'linear'
                    nSamples = p.Results.NSamples;
                    
                    % Normalized time vector
                    t = linspace(0, 1, nSamples)';
        
                    x = startPos(1) + (endPos(1) - startPos(1)) * t;
                    y = startPos(2) + (endPos(2) - startPos(2)) * t;
                    z = startPos(3) + (endPos(3) - startPos(3)) * t;
        
                    trajectory = [x, y, z];
                    
                case 'circle' 
                    S = p.Results.Params; 
                    radius = S.radius; 
                    angle = S.angle; 
                    direction = S.direction; 
                    
                    nSamples = p.Results.NSamples;
                    theta = linspace(0, deg2rad(angle), nSamples)';
                    
                    % Rotation direction
                    if strcmpi(direction, 'right')
                        theta = -theta;
                    end
        
                    % Circle center (based on startPos)
                    cx = startPos(1);
                    cy = startPos(2);
                    cz = startPos(3);
        
                    % Generate trajectory
                    x = cx + radius * cos(theta);
                    y = cy + radius * sin(theta);
                    z = cz * ones(size(theta));
        
                    trajectory = [x, y, z];
                    
                case 'curve'
                    waypoints = p.Results.waypoints; 
                    if isempty(waypoints)
                        error('Waypoints parameters must be defined to use curve mode.');
                    end
                    
                    nSamples = p.Results.NSamples;
        
                    % Include start, waypoints, and end positions
                    pts = [startPos; waypoints; endPos];
        
                    % Accumulated distance for correct parameterization
                    d = [0; cumsum(sqrt(sum(diff(pts).^2, 2)))];
        
                    % Sampling vector
                    t = linspace(0, d(end), nSamples);
        
                    % Spline interpolation
                    x = spline(d, pts(:,1), t)';
                    y = spline(d, pts(:,2), t)';
                    z = spline(d, pts(:,3), t)';
        
                    trajectory = [x, y, z];
            end  
        end
        
        function [thetas_est, phis_est, K_est] = estimar_DOAs_MUSIC( ...
            obj, M, N, d, fc, c, ...
            thetas_true, phis_true, ...
            N_snapshots, snr_db, ...
            coherence, ...
            thetas_grid, phis_grid, ...
            seed, estimador)

            % Description: Convenience method integrating UPA creation, signal 
            %              generation, and MUSIC DOA estimation.
            %
            % Inputs:
            %   M, N        - Array dimensions (number of elements in x and y)
            %   d           - Inter-element spacing
            %   fc, c       - Carrier frequency (Hz) and speed of light (m/s)
            %   thetas_true - True zenith angles of sources
            %   phis_true   - True azimuth angles of sources
            %   N_snapshots - Number of temporal snapshots
            %   snr_db      - Signal-to-Noise Ratio in dB
            %   coherence   - Correlation coefficient between sources
            %   thetas_grid - Search grid for zenith angles (radians)
            %   phis_grid   - Search grid for azimuth angles (radians)
            %   seed        - RNG seed for reproducibility
            %   estimador   - Instantiated DOA estimator object
            %
            % Outputs:
            %   thetas_est  - Estimated zenith angles
            %   phis_est    - Estimated azimuth angles
            %   K_est       - Estimated number of sources
            
            % Create array object
            array_vca = UniformPlanarArray(M, N, d, d, fc, c);
        
            % Initialize signal model
            signal_model = SignalModel(array_vca);
        
            % Generate signals
            [X, ~, ~] = signal_model.generateSignals( ...
                thetas_true, phis_true, ...
                N_snapshots, snr_db, ...
                'complex_sinusoid', 0.0, 'white', seed);
            
            % Estimate number of sources (optional step based on estimator)
            K_est = estimador.calculateNumSources(X, coherence, false);
            
            % Execute MUSIC algorithm
            [thetas_est, phis_est] = estimador.estimate(X, K_est, 0.05*pi/180, coherence);
        end
        
        function [x_next, P_next] = step_plkf_3d(obj, x_prev, P_prev, measurements, s_curr, T, process_q, sigma_angles)
            
            % Description: Single step of the standard Pseudo-Linear Kalman Filter (PLKF).
            %
            % Prediction:
            % $$\mathbf{x}_{k|k-1} = \mathbf{F}\mathbf{x}_{k-1|k-1}$$
            % $$\mathbf{P}_{k|k-1} = \mathbf{F}\mathbf{P}_{k-1|k-1}\mathbf{F}^T + \mathbf{Q}$$
            %
            % Inputs:
            %   x_prev       - Previous state vector [6x1]
            %   P_prev       - Previous covariance matrix [6x6]
            %   measurements - Current angles [theta; phi] (radians)
            %   s_curr       - Current observer position [sx; sy; sz]
            %   T            - Sampling period
            %   process_q    - Process noise covariance matrix [6x6]
            %   sigma_angles - Measurement variance [sigma_theta; sigma_phi]
            %
            % Outputs:
            %   x_next, P_next - Updated state and covariance
            
            %% Transition Matrix 
            F = [1 0 0 T 0 0;
                 0 1 0 0 T 0;
                 0 0 1 0 0 T;
                 0 0 0 1 0 0;
                 0 0 0 0 1 0;
                 0 0 0 0 0 1];
        
            %% Prediction Step
            x_pred = F * x_prev;
            P_pred = F * P_prev * F' + process_q;
        
            %% Pseudo-Linear Construction
            sx = s_curr(1); 
            sy = s_curr(2); 
            sz = s_curr(3);
            
            theta_k = measurements(1); 
            phi_k   = measurements(2);
        
            % Estimated distance and elevation for noise weighting
            rel_pred = x_pred(1:3) - s_curr;
            d_hat = norm(rel_pred);
            d_hat = max(d_hat, 1); % Protect against division by zero
            
            % Predicted elevation to correct azimuth noise
            phi_hat = asin(rel_pred(3) / d_hat);
        
            % Azimuth (Theta) components
            H_th = [sin(theta_k), -cos(theta_k), 0, 0, 0, 0];
            z_th = sx * sin(theta_k) - sy * cos(theta_k);
            R_th = (d_hat * cos(phi_hat))^2 * sigma_angles(1);
        
            % Elevation (Phi) components
            H_ph = [-sin(phi_k)*cos(theta_k), -sin(phi_k)*sin(theta_k), cos(phi_k), 0, 0, 0];
            z_ph = -sx * sin(phi_k) * cos(theta_k) - sy * sin(phi_k) * sin(theta_k) + sz * cos(phi_k);
            R_ph = (d_hat)^2 * sigma_angles(2);
        
            % Stack matrices
            H_k = [H_th; H_ph];
            z_k = [z_th; z_ph];
            R_k = diag([R_th, R_ph]);
        
            %% Update Step
            S_k = H_k * P_pred * H_k' + R_k;
            K_k = P_pred * H_k' / S_k;
            
            % Pseudo-linear innovation
            innovation = z_k - H_k * x_pred;
            
            x_next = x_pred + K_k * innovation;
            P_next = (eye(6) - K_k * H_k) * P_pred;
        end
        
        function [x_next, P_next] = step_uplkf_3d(obj, x_prev, P_prev, measurements, s_curr, F, process_q, sigma_angles)

            % Description: Single step of the Weighted Unbiased Pseudo-Linear Kalman Filter (WUPLKF).
            % Introduces a bias correction term in the measurement noise covariance
            % to mitigate the correlation between the pseudo-linear matrix H and the noise.
            %
            % Inputs:
            %   x_prev       - Previous state vector [6x1]
            %   P_prev       - Previous covariance matrix [6x6]
            %   measurements - Current angles [theta; phi] (radians)
            %   s_curr       - Current observer position [sx; sy; sz]
            %   F            - State transition matrix [6x6]
            %   process_q    - Process noise covariance matrix [6x6]
            %   sigma_angles - Measurement variance [sigma_theta; sigma_phi]
            %
            % Outputs:
            %   x_next       - Updated state vector [6x1]
            %   P_next       - Updated covariance matrix [6x6]

            
            %% Prediction Step
            x_pred = F * x_prev;
            P_pred = F * P_prev * F' + process_q;
        
            %% Geometry & Pseudo-Linear Construction ---
            sx = s_curr(1); 
            sy = s_curr(2); 
            sz = s_curr(3);
        
            theta_k = measurements(1); 
            phi_k   = measurements(2);
        
            rel_pred = x_pred(1:3) - s_curr;
            d_hat = max(norm(rel_pred), 1);
            phi_hat = asin(rel_pred(3) / d_hat);
        
            % Azimuth components
            H_th = [sin(theta_k), -cos(theta_k), 0, 0, 0, 0];
            z_th = sx * sin(theta_k) - sy * cos(theta_k);
        
            % Elevation components
            H_ph = [-sin(phi_k)*cos(theta_k), -sin(phi_k)*sin(theta_k), cos(phi_k), 0, 0, 0];
            z_ph = -sx * sin(phi_k) * cos(theta_k) - sy * sin(phi_k) * sin(theta_k) + sz * cos(phi_k);
        
            H_k = [H_th; H_ph];
            z_k = [z_th; z_ph];
        
            % Base projected noise
            R_th = (d_hat * cos(phi_hat))^2 * sigma_angles(1);
            R_ph = (d_hat)^2 * sigma_angles(2);
            R_k = diag([R_th, R_ph]);
        
            %%  UPLKF Bias Correction
            % Derivates with respect to the angles
            J_full_theta = [ cos(theta_k),            sin(theta_k),           0, 0, 0, 0;
                             sin(phi_k)*sin(theta_k), -sin(phi_k)*cos(theta_k), 0, 0, 0, 0 ];
            
            J_full_phi   = [ 0,                       0,                      0,           0, 0, 0;
                             -cos(phi_k)*cos(theta_k), -cos(phi_k)*sin(theta_k), -sin(phi_k), 0, 0, 0 ];
            
            % Compensation term
            C_k = sigma_angles(1) * (J_full_theta * P_pred * J_full_theta') + ...
                  sigma_angles(2) * (J_full_phi   * P_pred * J_full_phi');
            
            R_unbiased = R_k + C_k;
        
            %% Update Step
            S_k = H_k * P_pred * H_k' + R_unbiased;
            K_k = P_pred * H_k' / S_k;
        
            innovation = z_k - H_k * x_pred;
        
            x_next = x_pred + K_k * innovation;
            P_next = (eye(6) - K_k * H_k) * P_pred;
        end
        
        function [x_k, P_k] = ukf_step(obj, f_func, h_func, x_prev, P_prev, z_k, Q, R, params)

            % Description: Single step of the Unscented Kalman Filter (UKF).
            %
            % Inputs:
            %   f_func - Function handle for state transition f(x)
            %   h_func - Function handle for measurement model h(x)
            %   x_prev - Previous state estimate [n x 1]
            %   P_prev - Previous covariance estimate [n x n]
            %   z_k    - Current measurement [m x 1]
            %   Q      - Process noise covariance
            %   R      - Measurement noise covariance
            %   params - Struct with UKF tuning parameters (.n, .alpha, .kappa, .beta)
            %
            % Outputs:
            %   x_k, P_k - Updated state and covariance estimates

            
            %% Initialization & Parameters
            n = params.n;
            alpha = params.alpha;
            kappa = params.kappa;
            beta = params.beta; 
            
            % Scaling factor lambda
            lambda = alpha^2 * (n + kappa) - n;
            
            %% Weight Calculation
            Wm = zeros(1, 2*n + 1);
            Wc = zeros(1, 2*n + 1);
            
            % Weights for the central point (i = 0)
            Wm(1) = lambda / (n + lambda);
            Wc(1) = (lambda / (n + lambda)) + (1 - alpha^2 + beta);
            
            % Weights for the remaining points (i = 1 to 2n)
            scalar_weight = 1 / (2 * (n + lambda));
            Wm(2:end) = scalar_weight;
            Wc(2:end) = scalar_weight;
        
            %% Sigma Points Generation
            c = n + lambda;
            try
                S = chol(c * P_prev, 'lower'); 
            catch
                % Fallback for numerical instability (ensure positive definite)
                S = chol(c * P_prev + 1e-9 * eye(n), 'lower'); 
            end
            
            sigma_points = zeros(n, 2*n + 1);
            
            % Central point
            sigma_points(:, 1) = x_prev;
            
            % Displaced points
            for i = 1:n
                sigma_points(:, i+1)     = x_prev + S(:, i);
                sigma_points(:, n+i+1)   = x_prev - S(:, i);
            end
        
            %% Time Update (Prediction)
            sigma_points_pred = zeros(n, 2*n + 1);
            for i = 1 : (2*n + 1)
                sigma_points_pred(:, i) = f_func(sigma_points(:, i));
            end
            
            % Predicted state mean
            x_pred = zeros(n, 1);
            for i = 1 : (2*n + 1)
                x_pred = x_pred + Wm(i) * sigma_points_pred(:, i);
            end
            
            % Predicted state covariance
            P_pred = zeros(n, n);
            for i = 1 : (2*n + 1)
                diff_val = sigma_points_pred(:, i) - x_pred;
                P_pred = P_pred + Wc(i) * (diff_val * diff_val');
            end
            % Inject process noise
            P_pred = P_pred + Q; 
        
            %% Measurement Update
            m = length(z_k);
            
            % Propagate sigma points through measurement function
            Z_sigma = zeros(m, 2*n + 1);
            for i = 1 : (2*n + 1)
                Z_sigma(:, i) = h_func(sigma_points_pred(:, i)); 
            end
            
            % Predicted measurement mean
            z_pred = zeros(m, 1);
            for i = 1 : (2*n + 1)
                z_pred = z_pred + Wm(i) * Z_sigma(:, i);
            end
        
            %% Covariance & State Estimation
            % Innovation covariance
            P_zz = zeros(m, m);
            for i = 1 : (2*n + 1)
                diff_z = Z_sigma(:, i) - z_pred;
                P_zz = P_zz + Wc(i) * (diff_z * diff_z');
            end
            % Inject measurement noise
            P_zz = P_zz + R; 
            
            % Cross covariance
            P_xz = zeros(n, m);
            for i = 1 : (2*n + 1)
                diff_x = sigma_points_pred(:, i) - x_pred;
                diff_z = Z_sigma(:, i) - z_pred;
                P_xz = P_xz + Wc(i) * (diff_x * diff_z');
            end
            
            % Kalman Gain
            K = P_xz / P_zz; 
            
            % Final state update
            x_k = x_pred + K * (z_k - z_pred);
            
            % Final covariance update
            P_k = P_pred - K * P_zz * K';
        end
    end
end