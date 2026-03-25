classdef fast_MUSIC

    % Description: Fast MUSIC algorithm for 2D Direction of Arrival (DOA) 
    %              estimation utilizing a coarse-to-fine search grid approach 
    %              to reduce computational complexity.
    
    properties
        array                   % Array configuration object
        Step_Coarse = 2         % Coarse search step resolution (degrees)
        Step_Fine   = 0.05      % Fine search step resolution (degrees)
        Window_Fine = 5         % Fine search window span (+/- degrees)
        Theta_lim   = [0 90]    % Evaluation limits for elevation/zenith (degrees)
        Phi_lim     = [-180 180]% Evaluation limits for azimuth (degrees)
    end
    
    methods
        function obj = fast_MUSIC(array)

            % Description: Constructor for the fast_MUSIC class.
            %
            % Inputs:
            %   array - Antenna array configuration object
            %
            % Outputs:
            %   obj   - Initialized instance of the class

            obj.array = array;
        end
        
        function [thetas_est, phis_est, debug] = estimate(obj, X, K, coherent)

            % Description: Estimates the 2D DOAs applying a two-stage 
            %              (coarse and fine) MUSIC algorithm.
            %
            % Inputs:
            %   X          - Received signal matrix (Sensors x Snapshots)
            %   K          - Expected number of signal sources
            %   coherent   - Boolean indicating if spatial smoothing is required
            %
            % Outputs:
            %   thetas_est - Estimated elevation angles (radians)
            %   phis_est   - Estimated azimuth angles (radians)
            %   debug      - Struct containing intermediate grid and spectrum data

            M = obj.array.M;
            N = obj.array.N;
            Lx = floor(M/2);
            Ly = floor(N/2);
            
            U_n = obj.calculateNoiseSubspace(X, K, M, N, Lx, Ly, coherent);
            
            % Stage 1: Coarse Scan
            th_coarse = obj.Theta_lim(1) : obj.Step_Coarse : obj.Theta_lim(2);
            ph_coarse = obj.Phi_lim(1)   : obj.Step_Coarse : obj.Phi_lim(2);
            S_coarse = obj.calculateSpectrumGrid(th_coarse, ph_coarse, U_n, coherent);
            
            % Locate local maxima in the coarse spectrum
            BW = imregionalmax(S_coarse);
            [row, col] = find(BW);
            peak_vals = S_coarse(sub2ind(size(S_coarse), row, col));
            
            S_vec = S_coarse(:);
            median_S = median(S_vec);
            % sigma_S = std(S_vec); % (Not used in coarse thresholding, kept for reference if needed)
            
            % Define threshold (3 times the median)
            threshold = median_S * 3;
            keep_idx = peak_vals >= threshold;
            if any(~keep_idx)
                fprintf('Warning: Some coarse peaks are below %d dB and will be ignored.\n', abs(threshold));
            end
            
            row = row(keep_idx);
            col = col(keep_idx);
            peak_vals = peak_vals(keep_idx);
            
            % Sort peaks in descending order
            [~, idx_sort] = sort(peak_vals, 'descend');
            n_peaks = min(K, length(idx_sort));
            candidates_th = th_coarse(col(idx_sort(1:n_peaks)));
            candidates_ph = ph_coarse(row(idx_sort(1:n_peaks)));
            
            n_ignored = K - n_peaks;
            
            % Store stage 1 debug information
            debug.stage1.grid_th = th_coarse;
            debug.stage1.grid_ph = ph_coarse;
            debug.stage1.spectrum = S_coarse;
            debug.stage1.peaks = [candidates_th', candidates_ph'];
            
            % Stage 2: Fine Scan
            thetas_all = [];
            phis_all   = [];
            
            for k = 1:n_peaks
                center_th = candidates_th(k);
                center_ph = candidates_ph(k);
                
                th_fine = (center_th - obj.Window_Fine) : obj.Step_Fine : (center_th + obj.Window_Fine);
                ph_fine = (center_ph - obj.Window_Fine) : obj.Step_Fine : (center_ph + obj.Window_Fine);
                
                % Constrain fine grid to defined angular limits
                th_fine = th_fine(th_fine >= obj.Theta_lim(1) & th_fine <= obj.Theta_lim(2));
                ph_fine = ph_fine(ph_fine >= obj.Phi_lim(1) & ph_fine <= obj.Phi_lim(2));
                
                S_fine = obj.calculateSpectrumGrid(th_fine, ph_fine, U_n, coherent);
                
                % Locate local maxima in the fine spectrum
                BW_fine = imregionalmax(S_fine);
                [row_f, col_f] = find(BW_fine);
                peak_vals_f = S_fine(sub2ind(size(S_fine), row_f, col_f));
                
                median_S_fine = median(S_coarse(:));
                sigma_S_fine  = std(S_coarse(:));
                
                % Define threshold (median + 5 times standard deviation)
                threshold_fine = median_S_fine + 5 * sigma_S_fine;
    
                keep_idx_f = peak_vals_f >= threshold_fine;
                
                row_f = row_f(keep_idx_f);
                col_f = col_f(keep_idx_f);
                peak_vals_f = peak_vals_f(keep_idx_f);
                
                % Sort descending
                [~, idx_sort_f] = sort(peak_vals_f, 'descend');
                
                % Select the top candidates considering previously ignored peaks
                n_take = min(1 + n_ignored, length(idx_sort_f));
                r_final = row_f(idx_sort_f(1:n_take));
                c_final = col_f(idx_sort_f(1:n_take));
                
                % Convert indices back to physical angles (in radians)
                th_angles = deg2rad(th_fine(c_final));
                ph_angles = deg2rad(ph_fine(r_final));
                
                % Aggregate final peaks
                thetas_all = [thetas_all, th_angles];
                phis_all   = [phis_all, ph_angles];
                
                % Store stage 2 debug information
                debug.stage2(k).grid_th = th_fine;
                debug.stage2(k).grid_ph = ph_fine;
                debug.stage2(k).spectrum = S_fine;
                debug.stage2(k).peak = [rad2deg(th_angles)', rad2deg(ph_angles)'];
            end
            
            thetas_est = thetas_all;
            phis_est   = phis_all;
        end
        
        function spectrum = calculateSpectrumGrid(obj, theta_grid, phi_grid, U_n, coherent)

            % Description: Evaluates the MUSIC spectrum over a defined angular grid.
            %
            % Inputs:
            %   theta_grid - Evaluation grid for elevation angles (degrees)
            %   phi_grid   - Evaluation grid for azimuth angles (degrees)
            %   U_n        - Noise subspace eigenvector matrix
            %   coherent   - Boolean indicating spatial smoothing necessity
            %
            % Outputs:
            %   spectrum   - Computed 2D spectrum matrix

            [THETA_GRID, PHI_GRID] = meshgrid(theta_grid, phi_grid);
            theta_flat = THETA_GRID(:);
            phi_flat   = PHI_GRID(:);
            
            M = obj.array.M;
            N = obj.array.N;
            Lx = floor(M/2);
            Ly = floor(N/2);
            
            % Convert degrees to radians for array manifold calculation
            A = obj.array.arrayManifold(deg2rad(theta_flat).', deg2rad(phi_flat).'); 
            
            if coherent
                idx = obj.getSmoothingIndices(M, N, Lx, Ly);
                A = A(idx, :); 
            end
            
            projections = sum(abs(U_n' * A).^2, 1); 
            spectrum_flat = 1 ./ (projections + 1e-12);
            spectrum = reshape(spectrum_flat, length(phi_grid), length(theta_grid));
        end
        
        function K_est = calculateNumSources(obj, X, coherent, plt)

            % Description: Estimates the number of sources applying the Akaike 
            %              Information Criterion (AIC).
            %
            % Inputs:
            %   X        - Received signal matrix
            %   coherent - Boolean indicating spatial smoothing necessity
            %   plt      - Boolean to toggle plot generation
            %
            % Outputs:
            %   K_est    - Estimated number of signal sources

            arguments
                obj
                X
                coherent (1,1) logical = false 
                plt (1,1) logical = false 
            end
            
            % Step 1: Establish correct dimensions
            [~, n_snapshots] = size(X); 
            
            if coherent
                % Effective array size reduces when applying spatial smoothing
                M_rows = obj.array.M;
                N_cols = obj.array.N;
                Lx = floor(M_rows/1.75);
                Ly = floor(N_cols/1.75);
                R = fast_MUSIC.spatialSmoothing2D(X, M_rows, N_cols, Lx, Ly);
                
                % Virtual sensors count corresponds to the size of the smoothed matrix
                M_eff = size(R, 1); 
            else
                R = (X * X') / n_snapshots;
                M_eff = size(X, 1); 
            end
            
            % Step 2: Eigen decomposition
            [~, D] = eig(R);
            eigs_R = sort(abs(diag(D)), 'descend'); 
            
            % Step 3: Compute AIC sequence
            % Evaluation is limited to M_eff - 1
            aic_vals = zeros(M_eff - 1, 1);
            
            for k = 0:M_eff-2
                m = M_eff - k; 
                noise_eigs = eigs_R(k+1:end);
                
                % Geometric and arithmetic means calculated in logarithmic scale
                log_geom_mean = mean(log(noise_eigs + 1e-12)); 
                arith_mean = mean(noise_eigs);
                
                % Log-Likelihood corrected by snapshot count
                log_likelihood = m * n_snapshots * (log_geom_mean - log(arith_mean));
                
                % Penalty assignment for complex signals
                penalty = 2 * k * (2 * M_eff - k);
                
                aic_vals(k+1) = -2 * log_likelihood + penalty;
            end
            
            % Step 4: Final source estimation
            [~, min_idx] = min(aic_vals);
            K_est = min_idx - 1; 
            
            % Plot generation for debugging
            if plt
                figure; plot(0:M_eff-2, aic_vals, '-o');
                title('AIC Criterion by Number of Sources');
                xlabel('K (Sources)'); ylabel('AIC Value'); grid on;
            end
        end
        
        function spectrum_2d = musicSpectrumVectorized(obj, X, K, theta_grid, phi_grid, coherent)

            % Description: Computes the 2D MUSIC spectrum using matrix vectorization.
            %
            % Inputs:
            %   X          - Signal matrix
            %   K          - Number of sources
            %   theta_grid - Vector of evaluation angles for elevation
            %   phi_grid   - Vector of evaluation angles for azimuth
            %   coherent   - Boolean indicating spatial smoothing necessity
            %
            % Outputs:
            %   spectrum_2d - Computed 2D spectrum matrix

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
    end
    
    methods (Static, Access = private)
        function U_n = calculateNoiseSubspace(X, K, M, N, Lx, Ly, coherent)

            % Description: Extracts the noise subspace from the signal correlation matrix.
            %
            % Inputs:
            %   X        - Signal matrix
            %   K        - Number of target sources
            %   M, N     - Total array dimensions
            %   Lx, Ly   - Sub-array dimensions for smoothing
            %   coherent - Boolean determining if spatial smoothing is applied
            %
            % Outputs:
            %   U_n      - Noise subspace eigenvector matrix

            if coherent
                R = fast_MUSIC.spatialSmoothing2D(X, M, N, Lx, Ly);
            else
                R = (X * X') / size(X, 2);
            end
            [V, D] = eig(R);
            [~, idx] = sort(diag(D), 'descend');
            V = V(:, idx);
            U_n = V(:, K+1:end);
        end
        
        function Rss = spatialSmoothing2D(X, Nx, Ny, Lx, Ly)

            % Description: Performs 2D Spatial Smoothing to decorrelate coherent signals.
            %
            % Inputs:
            %   X      - Signal matrix
            %   Nx, Ny - Total planar array dimensions
            %   Lx, Ly - Sub-array dimensions
            %
            % Outputs:
            %   Rss    - Smoothed spatial covariance matrix

            Nant = Nx * Ny;
            snaps = size(X, 2);
            
            if size(X,1) ~= Nant
                error('Inconsistent dimension: Matrix X must have Nx*Ny rows.');
            end
            
            num_sub_x = Nx - Lx + 1;
            num_sub_y = Ny - Ly + 1;
            num_subarrays = num_sub_x * num_sub_y;
            
            Rss = zeros(Lx*Ly, Lx*Ly);
            
            for ix = 1:num_sub_x
                for iy = 1:num_sub_y
                    idx = [];
                    for x = ix:(ix+Lx-1)
                        base = (x-1)*Ny;
                        idx = [idx, base + iy : base + iy + Ly - 1];
                    end
                    Xsub = X(idx, :);
                    Rsub = (Xsub * Xsub') / snaps;
                    Rss = Rss + Rsub;
                end
            end
            
            Rss = Rss / num_subarrays;
        end
        
        function indices = getSmoothingIndices(Nx, Ny, Lx, Ly)

            % Description: Constructs the index mapping for 2D spatial smoothing limits.
            %
            % Inputs:
            %   Nx, Ny - Total dimensions of the main array
            %   Lx, Ly - Dimensions corresponding to the sub-array
            %
            % Outputs:
            %   indices - Numeric mapping of the sub-array selection

            indices = [];
            ix = 1;
            iy = 1;
            for x = ix:(ix+Lx-1)
                base = (x-1)*Ny;
                indices = [indices, base + iy : base + iy + Ly - 1];
            end
        end
    end
end