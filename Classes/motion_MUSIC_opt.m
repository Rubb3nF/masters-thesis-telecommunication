classdef motion_MUSIC_opt < handle
    properties
        array  
    end
    
    properties (Access = public)
        history = {}; 
        
        % Configuration
        search_radius_tracking = 10;
        
        % State Machine
        state = 0; % 0: Normal, 1: Collision, 2: Separation
        separation_threshold = 10*pi/180; 
        
        % Debug info
        last_pred_theta
        last_pred_phi
        last_r_theta
        last_r_phi
        
        fast_solver
    end
    
    methods
        function obj = motion_MUSIC_opt(array_input)

            % Constructor: Initializes the motion_MUSIC_opt object.
            %
            % Inputs:
            %   array_input - Antenna array configuration object
            %
            % Outputs:
            %   obj         - Initialized instance of the class

            obj.array = array_input;
            obj.history = {}; 
            obj.fast_solver = fast_MUSIC(array_input);
            obj.last_pred_theta = []; obj.last_pred_phi = [];
            obj.last_r_theta = [];    obj.last_r_phi = [];
        end
        
        function [thetas, phis] = estimate(obj, X, K, safety_guard_rad, coherent)

            % Description: Estimates the Direction of Arrival (DoA) of K sources with
            %              tracking capabilities and a state machine for collisions.
            %
            % Inputs:
            %   X                - Received signal matrix (Sensors x Snapshots)
            %   K                - Number of sources to track
            %   safety_guard_rad - Safety guard radius for kinematic tracking
            %   coherent         - Boolean indicating if spatial smoothing
            %                      is applied
            %
            % Outputs:
            %   thetas           - Estimated elevation angles (radians)
            %   phis             - Estimated azimuth angles (radians)

            if nargin < 4, error('Missing arguments'); end
        
            % Phase 1: Initialization
            % Verify if a valid history exists for K sources
            valid_cols = 0;
            if ~isempty(obj.history)
                valid_cols = sum(sum(~cellfun(@isempty, obj.history)) >= 3);
            end
            
            if valid_cols < K
                [thetas, phis] = obj.fast_solver.estimate(X, K, coherent);
                if ~isempty(thetas)
                    obj.clustering([thetas(:), phis(:)]); 
                end
                return; 
            end
        
            % Phase 2: Prediction instruments
            prev3 = vertcat(obj.history{1, :}); % t-3
            prev2 = vertcat(obj.history{2, :}); % t-2
            prev1 = vertcat(obj.history{3, :}); % t-1
            
            pred_t = zeros(K, 1); pred_p = zeros(K, 1);
            r_thetas = zeros(K, 1); r_phis = zeros(K, 1);
            
            for k = 1:K
                % Kinematic prediction
                [pred_t(k), pred_p(k)] = obj.predictNextPositionRodrigues(...
                    prev3(k,1), prev3(k,2), prev2(k,1), prev2(k,2), prev1(k,1), prev1(k,2));
                
                % Adaptive radius computation
                [r_thetas(k), r_phis(k)] = obj.calculateAdaptiveRadii(...
                    prev2(k,1), prev2(k,2), prev1(k,1), prev1(k,2),pred_t(k), pred_p(k), pred_t(k), safety_guard_rad);
            end
            
            obj.last_pred_theta = pred_t; obj.last_pred_phi = pred_p;
            obj.last_r_theta = r_thetas;  obj.last_r_phi = r_phis;
            
            % Phase 2.5: Global threshold estimation for coarse-fine search
            coarse_step = 5; 
            [theta_coarse, phi_coarse] = meshgrid(0:coarse_step:90, -180:coarse_step:179);
            
            coarse_spectrum = musicSpectrumPointList(obj, X, K, theta_coarse(:), phi_coarse(:), coherent);
            
            global_median = median(coarse_spectrum);
            global_sigma   = std(coarse_spectrum);
            
            global_threshold = global_median + 3 * global_sigma;
            
            % Phase 3: Measurement (Windowed search)
            detected_thetas = zeros(K, 1);
            detected_phis   = zeros(K, 1);
            peak_powers     = zeros(K, 1);
            found_flags     = false(K, 1);
            
            for k = 1:K
                [theta_est, phi_est, peak_val] = obj.findPeakOptimized(...
                X, K, pred_t(k), pred_p(k), r_thetas(k), r_phis(k), coherent);
                
                if peak_val > global_threshold
                    detected_thetas(k) = theta_est;
                    detected_phis(k)   = phi_est;
                    peak_powers(k)     = peak_val;
                    found_flags(k)     = true;
                else
                    % If the threshold is not exceeded, consider the source lost or overlapped
                    found_flags(k)     = false;
                    fprintf('Source %d: Detected peak (%.2f) is statistical noise (Threshold: %.2f)\n', ...
                            k, peak_val, global_threshold);
                end
            end
            
            % Phase 4: State machine
            % Convert angles to 3D unit vectors to avoid cuts at +/- 180
            U_det = [sin(detected_thetas).*cos(detected_phis), ...
                     sin(detected_thetas).*sin(detected_phis), ...
                     cos(detected_thetas)];
            
            % Round to 3 decimals (approx 0.05 degrees tolerance)
            U_det_round = round(U_det, 3);
            unique_peaks_3d = unique(U_det_round, 'rows');
            num_unique = size(unique_peaks_3d, 1);
            
            switch obj.state
                case 0 % NORMAL
                    if num_unique < K && all(found_flags)
                        obj.state = 1; 
                        fprintf('State change: NORMAL -> COLLISION\n');
                    end
                case 1 % COLLISION
                    if num_unique == K
                        obj.state = 2;
                        fprintf('State change: COLLISION -> SEPARATION\n');
                    end
                case 2 % SEPARATION
                    dists = pdist(unique_peaks_3d); % 3D Euclidean distance
                    
                    % Adjust threshold to chordal distance: d = 2*sin(theta/2)
                    chordal_threshold = 2 * sin(obj.separation_threshold / 2);
                    
                    if min(dists) > chordal_threshold
                        obj.state = 0; 
                        fprintf('State change: SEPARATION -> NORMAL\n');
                    elseif num_unique < K
                        obj.state = 1;
                        fprintf('State change: SEPARATION -> COLLISION\n');
                    end
            end
            
            % Catastrophic failure check (signal lost)
            if any(~found_flags) && obj.state ~= 1
                 fprintf("Tracking failed (Signal lost). Resetting...\n");
                 tic
                 [thetas, phis] = obj.fast_solver.estimate(X, K, coherent);
                 toc
                 % History reset
                 obj.history = {};
                 if ~isempty(thetas)
                     obj.clustering([thetas(:), phis(:)]); 
                     obj.state = 0;
                 end
                 return;
            end
            
            % Phase 5: Differentiated update
            if obj.state == 1
                % COLLISION MODE
                % Ignore measured parameters for history and utilize prediction inertia.
                % Prevents shared peak noise from corrupting velocity.
                thetas = pred_t; 
                phis   = pred_p;
                
                % Manual history update (bypass clustering)
                % Assumes prediction k corresponds to source k
                for k = 1:K
                    obj.history(1:end-1, k) = obj.history(2:end, k);
                    obj.history{3, k} = [pred_t(k), pred_p(k)];
                end
                
            else
                % NORMAL / SEPARATION MODE
                % Utilize measured parameters
                thetas = detected_thetas;
                phis   = detected_phis;
                
                % Apply clustering for correct assignment
                obj.clustering([thetas(:), phis(:)]);
            end
        end
        
        function clustering(obj, new_pairs)

            % Description: Assigns new angle pairs to historical trajectories based 
            %              on 3D Euclidean distances and coherence metrics.
            %
            % Inputs:
            %   new_pairs - Nx2 matrix containing [theta, phi] pairs of new detections
            %
            % Outputs:
            %   None (Updates the internal history buffer)

            if isempty(new_pairs), return; end
            MAX_ROWS = 3;
            [~, existing_K] = size(obj.history);
            num_new = size(new_pairs, 1);
    
            % Initialize if history is empty
            if existing_K == 0
                for i = 1:num_new, obj.history{1, i} = new_pairs(i, :); end
                return;
            end
            
            % Dynamic source growth handling
            % If there are more peaks than columns, generate new columns
            if num_new > existing_K
                for i = (existing_K+1):num_new
                    obj.history{1, i} = zeros(1,length(new_pairs(i, :)));
                end
                existing_K = num_new; 
            end
    
            % Step 1: Extract previous data in 3D coordinates
            U_last = zeros(existing_K, 3);
            U_vel = zeros(existing_K, 3);
            for k = 1:existing_K
                idx = find(~cellfun(@isempty, obj.history(:, k)), 1, 'last');
                if ~isempty(idx)
                    t_u = obj.history{idx, k}(1); p_u = obj.history{idx, k}(2);
                    U_last(k, :) = [sin(t_u)*cos(p_u), sin(t_u)*sin(p_u), cos(t_u)];
                    
                    if idx > 1
                         t_prev = obj.history{idx-1, k}(1); p_prev = obj.history{idx-1, k}(2);
                         U_prev = [sin(t_prev)*cos(p_prev), sin(t_prev)*sin(p_prev), cos(t_prev)];
                         U_vel(k, :) = U_last(k, :) - U_prev; % Velocity vector
                    end
                else
                    U_last(k, :) = [inf, inf, inf];
                end
            end
    
            % Step 2: Convert new angle pairs to 3D coordinates
            U_new = zeros(num_new, 3);
            for i = 1:num_new
                t_n = new_pairs(i,1); p_n = new_pairs(i,2);
                U_new(i,:) = [sin(t_n)*cos(p_n), sin(t_n)*sin(p_n), cos(t_n)];
            end
    
            % Step 3: Compute cost matrix strictly in 3D
            cost_mat = zeros(num_new, existing_K);
            for i = 1:num_new
                for k = 1:existing_K
                    vec_prop = U_new(i,:) - U_last(k,:);
                    dist_euc = norm(vec_prop); % Authentic chordal distance
                    
                    if (obj.state == 2 || obj.state == 0) && norm(U_vel(k,:)) > 1e-4
                        norm_p = norm(vec_prop);
                        if norm_p > 1e-4
                            cohere = dot(U_vel(k,:), vec_prop) / (norm(U_vel(k,:)) * norm_p);
                            factor = 2 - cohere; 
                            cost_mat(i,k) = dist_euc * factor;
                        else
                            cost_mat(i,k) = dist_euc;
                        end
                    else
                        cost_mat(i,k) = dist_euc;
                    end
                end
            end
            
           [~, assignments_k] = min(cost_mat, [], 2);
           
            % Step 4: Update the history buffer
            for i = 1:num_new
                k_dest = assignments_k(i);
                current_pair = new_pairs(i, :);
                col_data = obj.history(:, k_dest);
                n_occupied = sum(~cellfun(@isempty, col_data));
                
                if n_occupied < MAX_ROWS
                    obj.history{n_occupied + 1, k_dest} = current_pair;
                else
                    obj.history(1:end-1, k_dest) = obj.history(2:end, k_dest);
                    obj.history{MAX_ROWS, k_dest} = current_pair;
                end
            end
        end
        
        function [theta_est, phi_est, peak_val] = findPeakOptimized(obj, X, K, t_center, p_center, r_theta, r_phi, coherent)

            % Description: Performs a localized numerical optimization to find the 
            %              peak MUSIC spectrum value within an elliptical region.
            %
            % Inputs:
            %   X        - Signal matrix
            %   K        - Number of sources
            %   t_center - Elevation center of the search area (radians)
            %   p_center - Azimuth center of the search area (radians)
            %   r_theta  - Search radius in elevation
            %   r_phi    - Search radius in azimuth
            %   coherent - Boolean indicating spatial smoothing necessity
            %
            % Outputs:
            %   theta_est - Estimated elevation
            %   phi_est   - Estimated azimuth
            %   peak_val  - MUSIC spectrum magnitude at the peak

            M_val = obj.array.M; 
            N_val = obj.array.N; 
            
            % Step 1: Subspace preparation
            if coherent
                % When coherent spatial smoothing is desired, assuming regular UPA
                Lx = floor(M_val/2); Ly = floor(N_val/2);
                U_n = obj.calculateNoiseSubspace(X, K, M_val, N_val, Lx, Ly, coherent);
                
                % Virtual grid generation for the smoothed subarray
                vec_x = (0:Lx-1) - (Lx-1)/2;
                vec_y = (0:Ly-1) - (Ly-1)/2;
                d_spacing = obj.array.dx; 
                [GridX, GridY] = ndgrid(vec_x, vec_y);
                px = GridX(:) * d_spacing; 
                py = GridY(:) * d_spacing;
                
            else
                % General mode: random array or full UPA
                % Calculate subspace using all available antennas
                U_n = obj.calculateNoiseSubspace(X, K, M_val, N_val, M_val, N_val, false);
                
                % Obtain physical positions from the array configuration object
                pos_meters = obj.array.element_positions;
                lambda = obj.array.wavelength;
                
                % Extract coordinates and normalize by wavelength
                px = (pos_meters(1, :).') / lambda; 
                py = (pos_meters(2, :).') / lambda;
            end
            
            % Dimensional integrity verification
            if length(px) ~= size(U_n, 1)
                error('Dimension error: The number of positions (%d) does not match the subspace (%d). If it is a Random array, ensure coherent is set to false.', length(px), size(U_n, 1));
            end
        
            % Noise projection matrix
            P_noise = U_n * U_n'; 
            
            % Step 2: Configure optimizer
            cost_func = @(ang) obj.internal_music_cost_oval(ang, px, py, P_noise, t_center, p_center, r_theta, r_phi);
        
            start_point = [t_center, p_center];
            options = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6);
            
            % Step 3: Execute optimization (fminsearch)
            [final_pos, min_cost] = fminsearch(cost_func, start_point, options);
            
            theta_est = final_pos(1);
            phi_est   = final_pos(2);
            
            % Critical normalization to ensure azimuth remains strictly within bounds
            phi_est = mod(phi_est + pi, 2*pi) - pi;
            
            peak_val = 1 / (min_cost + eps); 
        end
        
        function val = internal_music_cost_oval(obj, ang, px, py, P_noise, tc, pc, rt, rp)

            % Description: Internal cost function evaluating the MUSIC spectrum 
            %              inversely, restricted within an elliptical boundary.
            %
            % Inputs:
            %   ang     - [theta, phi] coordinates to evaluate
            %   px, py  - Antenna array coordinate vectors
            %   P_noise - Noise subspace projection matrix
            %   tc, pc  - Center coordinates of the ellipse
            %   rt, rp  - Radii of the ellipse
            %
            % Outputs:
            %   val     - Cost value (high values outside boundaries)

            th = ang(1);
            ph = ang(2);
            
            % Elliptical restriction evaluation
            % Calculate azimuth difference and constrain it within the [-pi, pi] range
            diff_ph = ph - pc;
            diff_ph = mod(diff_ph + pi, 2*pi) - pi; 
            
            % Normalized distance to the center (ellipse equation)
            dist_norm = ((th - tc) / rt)^2 + (diff_ph / rp)^2;
            
            if dist_norm > 1.0
                val = 1e10 * dist_norm; 
                return;
            end
            
            % Physical domain restriction
            if th < 0 || th > pi/2
                val = 1e10; 
                return;
            end
        
            % Standard MUSIC projection calculation
            u = sin(th) * cos(ph);
            v = sin(th) * sin(ph);
            
            phase = 2 * pi * (px * u + py * v);
            a = exp(1j * phase);
            
            val = real(a' * P_noise * a);
        end
        
        function [theta_valid, phi_valid] = generateMultipleSearchGrids(obj, pred_t, pred_p, r_thetas, r_phis, grid_step)

            % Description: Generates an aggregated set of valid grid points enclosed 
            %              by multiple search ellipses.
            %
            % Inputs:
            %   pred_t    - Kx1 vector of elevation centers
            %   pred_p    - Kx1 vector of azimuth centers
            %   r_thetas  - Kx1 vector of elevation radii
            %   r_phis    - Kx1 vector of azimuth radii
            %   grid_step - Spacing resolution for the grid
            %
            % Outputs:
            %   theta_valid - Column vector of valid evaluation points (elevation)
            %   phi_valid   - Column vector of valid evaluation points (azimuth)

            K = length(pred_t);
            all_theta = [];
            all_phi = [];
            
            for k = 1:K
                tc = pred_t(k);
                pc = pred_p(k);
                rt = r_thetas(k);
                rp = r_phis(k);
                
                % Step 1: Define boundary limits enveloping oval k
                t_min = max(tc - rt, 0); 
                t_max = min(tc + rt, pi/2);
                
                % A minor margin is implicitly handled to prevent edge loss due to rounding
                theta_vec = t_min : grid_step : t_max;
                phi_vec   = (pc - rp) : grid_step : (pc + rp);
                
                if isempty(theta_vec) || isempty(phi_vec), continue; end
                
                % Step 2: Construct localized mesh
                [TG, PG] = meshgrid(theta_vec, phi_vec);
                
                % Step 3: Apply elliptical mask
                % Distance is normalized by the respective adaptive radii
                term_theta = ((TG - tc) / rt).^2;
                term_phi   = ((PG - pc) / rp).^2;
                mask = (term_theta + term_phi) <= 1.0;
                
                % Step 4: Aggregate enclosed coordinates
                all_theta = [all_theta; TG(mask)];
                all_phi   = [all_phi; PG(mask)];
            end
            
            % Step 5: Consolidate data by removing duplicate coordinates
            % A slight rounding is applied to ensure uniqueness filters out numerical noise
            points = unique(round([all_theta, all_phi] / (grid_step/10)) * (grid_step/10), 'rows');
            
            theta_valid = points(:, 1);
            phi_valid   = points(:, 2);
        end
        
        function spectrum_2d = musicSpectrumVectorized(obj, X, K, theta_grid, phi_grid, coherent)

            % Description: Computes a comprehensive 2D MUSIC spectrum across a predefined grid.
            %
            % Inputs:
            %   X          - Signal matrix
            %   K          - Number of sources
            %   theta_grid - Vector of evaluation angles for elevation
            %   phi_grid   - Vector of evaluation angles for azimuth
            %   coherent   - Boolean indicating spatial smoothing necessity
            %
            % Outputs:
            %   spectrum_2d - Computed spectrum mapped to the dimensions of the grid

            M = obj.array.M;
            N = obj.array.N;
            Lx = floor(M/2);
            Ly = floor(N/2);
            U_n = obj.calculateNoiseSubspace(X, K, M, N, Lx, Ly, coherent);
            
            [THETA_GRID, PHI_GRID] = meshgrid(theta_grid, phi_grid);
            
            theta_flat = THETA_GRID(:);
            phi_flat = PHI_GRID(:);
            
            A = obj.array.arrayManifold(theta_flat.', phi_flat.'); 
            
            if coherent
                idx = obj.getSmoothingIndices(M, N, Lx, Ly);
                A = A(idx, :); 
            end
            
            projections = sum(abs(A' * U_n).^2, 2);
            
            spectrum_flat = 1 ./ (projections + 1e-12);
            
            spectrum_2d = reshape(spectrum_flat, length(phi_grid), length(theta_grid));
        end
        
       function spectrum_list = musicSpectrumPointList(obj, X, K, theta_valid, phi_valid, coherent)

            % Description: Computes the MUSIC spectrum exclusively at specific sparse points.
            %
            % Inputs:
            %   X           - Signal matrix
            %   K           - Number of sources
            %   theta_valid - Nx1 vector containing elevation test points
            %   phi_valid   - Nx1 vector containing azimuth test points
            %   coherent    - Boolean indicating spatial smoothing necessity
            %
            % Outputs:
            %   spectrum_list - Nx1 vector of evaluated spectrum magnitudes

            M = obj.array.M;
            N = obj.array.N;
            Lx = floor(M/2);
            Ly = floor(N/2);
            
            % Step 1: Subspace extraction
            U_n = obj.calculateNoiseSubspace(X, K, M, N, Lx, Ly, coherent);
            
            % Step 2: Manifold array generation directly from row vectors
            A = obj.array.arrayManifold(theta_valid.', phi_valid.'); 
            
            if coherent
                idx = obj.getSmoothingIndices(M, N, Lx, Ly);
                A = A(idx, :); 
            end
            
            % Step 3: Noise subspace projection
            projections = sum(abs(A' * U_n).^2, 2);
            
            % Step 4: Spectrum aggregation into a linear array structure
            spectrum_list = 1 ./ (projections + 1e-12);
        end
        
        function [theta_next, phi_next] = predictNextPositionRodrigues( ...
            obj, theta_prev2, phi_prev2, theta_prev1, phi_prev1, theta_curr, phi_curr)

            % Description: Estimates the upcoming spatial coordinate leveraging Rodrigues'
            %              rotation formula and historical kinematic trends.
            %
            % Inputs:
            %   theta_prev2, phi_prev2 - Angle states at t-2
            %   theta_prev1, phi_prev1 - Angle states at t-1
            %   theta_curr, phi_curr   - Current angle states at t
            %
            % Outputs:
            %   theta_next - Predicted elevation angle
            %   phi_next   - Predicted azimuth angle

            beta = 0.7;                 % Acceleration factor (0 = none, 1 = maximum)
            max_step = 45*pi/180;       % Maximum allowable angular step
            eps_axis = 1e-8;
        
            t0 = theta_prev2; p0 = phi_prev2;
            t1 = theta_prev1; p1 = phi_prev1;
            t2 = theta_curr;  p2 = phi_curr;
        
            v0 = [sin(t0)*cos(p0); sin(t0)*sin(p0); cos(t0)];
            v1 = [sin(t1)*cos(p1); sin(t1)*sin(p1); cos(t1)];
            v2 = [sin(t2)*cos(p2); sin(t2)*sin(p2); cos(t2)];
        
            % Compute initial rotational axis (0 to 1)
            k1 = cross(v0, v1);
            n1 = norm(k1);
        
            if n1 < eps_axis
                omega1 = [0;0;0];
            else
                k1 = k1 / n1;
                a1 = acos(max(min(dot(v0,v1),1),-1));
                omega1 = a1 * k1;
            end
        
            % Compute subsequent rotational axis (1 to 2)
            k2 = cross(v1, v2);
            n2 = norm(k2);
        
            if n2 < eps_axis
                omega2 = [0;0;0];
            else
                k2 = k2 / n2;
                a2 = acos(max(min(dot(v1,v2),1),-1));
                omega2 = a2 * k2;
            end
        
            % Perform angular velocity forecasting
            omega3 = omega2 + beta*(omega2 - omega1);
            alpha3 = norm(omega3);
        
            % Enforce security constraints
            alpha3 = min(alpha3, max_step);
        
            if alpha3 < eps_axis
                v3 = v2;
            else
                k3 = omega3 / norm(omega3);
        
                % Complete Rodrigues rotation
                v3 = v2*cos(alpha3) + ...
                     cross(k3,v2)*sin(alpha3) + ...
                     k3*dot(k3,v2)*(1-cos(alpha3));
        
                v3 = v3 / norm(v3);  % Defensive structural normalization
            end
        
            % Translate resultant vector back to spherical coordinates
            theta_next = acos(v3(3));
            phi_next   = atan2(v3(2), v3(1));
        end
        
        function [r_theta, r_phi] = calculateAdaptiveRadii(obj, theta_prev2, phi_prev2, theta_prev, phi_prev, theta_curr, phi_curr, theta_curr_rad, grid_step)
            % Description: Computes dynamic search radii proportional to spatial velocity.
            %
            % Inputs:
            %   theta_prev2, phi_prev2 - Angle states at t-2
            %   theta_prev, phi_prev   - Previous angle state
            %   theta_curr, phi_curr   - Current angle state
            %   theta_curr_rad         - Current reference elevation
            %   grid_step              - Precision grid step configuration
            %
            % Outputs:
            %   r_theta - Computed search radius for elevation
            %   r_phi   - Computed search radius for azimuth
            K_safety = 1;      % Velocity multiplier margin factor
            beta = 0.7;        % Acceleration weighting factor
            R_min_base = 1*pi/180;    % Minimum allowed angular expansion
            
            % Compute baseline absolute variation in elevation
            v_theta = abs(theta_curr - theta_prev);
            v_theta_prev = abs(theta_prev - theta_prev2);
            a_theta = abs(v_theta - v_theta_prev);
            
            % Process variation in azimuth, counteracting discontinuous jumps across +/- 180
            diff_phi = phi_curr - phi_prev;
            v_phi = abs(mod(diff_phi + pi, 2*pi) - pi); 
            diff_phi_prev = phi_prev - phi_prev2;
            v_phi_prev = abs(mod(diff_phi_prev + pi, 2*pi) - pi);
            a_phi = abs(v_phi - v_phi_prev);
        
            % Apply expansions
            r_theta = (v_theta + beta * a_theta) * K_safety;
            if r_theta < (R_min_base + grid_step)
                r_theta = R_min_base + grid_step;
            end
            
            % Extract geometrical compensation metric based on elevation
            % Adjusts naturally large physical coverage parameters near polar extremes
            geom_factor = 1 / max(sin(theta_curr_rad), 0.1); 
            
            r_phi = (v_phi + beta * a_phi) * K_safety;
            
            if r_phi < ((R_min_base * geom_factor) + grid_step)
                r_phi = (R_min_base * geom_factor) + grid_step;
            end
            
            % Bound values against configured physical and computational maximums
            r_theta = min(max(r_theta, 1*pi/180), 10*pi/180);   
            r_phi   = min(max(r_phi, 1*pi/180), 20*pi/180);    
        end
        
        function [K_est, R] = calculateNumSources(obj, X, coherent, plt)

            % Description: Estimates the number of sources implicitly via the Akaike 
            %              Information Criterion.
            %
            % Inputs:
            %   X        - Received signal matrix
            %   coherent - Logical configuration flag for spatial smoothing usage
            %   plt      - Logical configuration flag to visualize output
            %
            % Outputs:
            %   K_est    - Estimated number of signal sources
            %   R        - Signal covariance matrix applied during estimation

            arguments
                obj
                X
                coherent (1,1) logical = false 
                plt (1,1) logical = false 
            end
            
            % Step 1: Matrix and dimensions definition
            [~, n_snapshots] = size(X); 
            
            if coherent
                M_rows = obj.array.M;
                N_cols = obj.array.N;
                Lx = floor(M_rows/1.75);
                Ly = floor(N_cols/1.75);
                R = motion_MUSIC_opt.spatialSmoothing2D(X, M_rows, N_cols, Lx, Ly);
                
                M_eff = size(R, 1); 
            else
                R = (X * X') / n_snapshots;
                M_eff = size(X, 1); 
            end
            
            % Step 2: Eigen decomposition
            [~, D] = eig(R);
            eigs_R = sort(abs(diag(D)), 'descend'); 
            
            % Step 3: Information criterion iterations
            aic_vals = zeros(M_eff - 1, 1);
            
            for k = 0:M_eff-2
                m = M_eff - k; 
                noise_eigs = eigs_R(k+1:end);
                
                % Numerical stability is guaranteed through geometric means with logarithmic approaches
                log_geom_mean = mean(log(noise_eigs + 1e-12)); 
                arith_mean = mean(noise_eigs);
                
                log_likelihood = m * n_snapshots * (log_geom_mean - log(arith_mean));
                penalty = 2 * k * (2 * M_eff - k);
                
                aic_vals(k+1) = -2 * log_likelihood + penalty;
            end
            
            % Step 4: Minimum assessment
            [~, min_idx] = min(aic_vals);
            K_est = min_idx - 1; 
            
            % Execute debug visualization upon request
            if plt
                figure; plot(0:M_eff-2, aic_vals, '-o');
                title('Akaike Information Criterion (AIC)');
                xlabel('K (Sources)'); ylabel('AIC Value'); grid on;
            end
        end
    end
    
    methods (Static, Access = private)
        function U_n = calculateNoiseSubspace(X, K, Nx, Ny, Lx, Ly, coherent)

            % Description: Extracts the noise subspace representation given K signals.
            %
            % Inputs:
            %   X        - Signal matrix
            %   K        - Estimated target count
            %   Nx, Ny   - Total array structural dimensions
            %   Lx, Ly   - Smooth spatial sub-array dimensions
            %   coherent - Boolean requirement switch
            %
            % Outputs:
            %   U_n      - Matrix encompassing noise eigenvectors

            if coherent
                R = motion_MUSIC_opt.spatialSmoothing2D(X, Nx, Ny, Lx, Ly);
            else
                R = (X * X') / size(X, 2);
            end
            [V, D] = eig(R);
            [~, idx] = sort(diag(D), 'descend');
            V = V(:, idx);
            
            U_n = V(:, K+1:end);
        end
        
        function Rss = spatialSmoothing2D(X, Nx, Ny, Lx, Ly)

            % Description: Decorrelates sources through sub-array processing.
            %
            % Inputs:
            %   X      - Target raw signal
            %   Nx, Ny - Native element dimensions
            %   Lx, Ly - Sub-array element dimensions
            %
            % Outputs:
            %   Rss    - Reconstructed spatially averaged covariance matrix

            Nant = Nx * Ny;
            snaps = size(X, 2);
        
            % Data validation procedures
            if size(X,1) ~= Nant
                error('Dimension inconsistency: Matrix X must correlate with total Nx*Ny rows.');
            end
        
            num_sub_x = Nx - Lx + 1;
            num_sub_y = Ny - Ly + 1;
            num_subarrays = num_sub_x * num_sub_y;
        
            Rss = zeros(Lx * Ly, Lx * Ly);
        
            % Perform matrix sub-sampling sweeps
            for ix = 1:num_sub_x
                for iy = 1:num_sub_y
                    
                    idx = [];
                    for x = ix:(ix + Lx - 1)
                        base = (x-1) * Ny;
                        idx = [idx, base + iy : base + iy + Ly - 1];
                    end
        
                    Xsub = X(idx, :);                  
                    Rsub = (Xsub * Xsub') / snaps;     
                    
                    Rss = Rss + Rsub;
                end
            end
        
            % Output mean consolidation
            Rss = Rss / num_subarrays;
        end
        
        function indices = getSmoothingIndices(Nx, Ny, Lx, Ly)

            % Description: Resolves relative index locations for reference overlapping subsets.
            %
            % Inputs:
            %   Nx, Ny - Total domain structure dimensions
            %   Lx, Ly - Reference subset structure sizes
            %
            % Outputs:
            %   indices - Numeric structural coordinates mapping mapping

            indices = [];
            ix = 1; % X-axis origin point
            iy = 1; % Y-axis origin point
            for x = ix:(ix+Lx-1)
                base = (x-1) * Ny;
                indices = [indices, base + iy : base + iy + Ly - 1];
            end
        end
    end
end