classdef MUSIC

    % Description: MUSIC algorithm for 1D or 2D Direction of Arrival (DOA) 
    %              estimation.
   
    properties
        array  % Array configuration object (e.g., UniformPlanarArray)
    end
    
    methods
        function obj = MUSIC(array)

            % Description: Constructor for the MUSIC class.
            %
            % Inputs:
            %   array - Antenna array configuration object
            %
            % Outputs:
            %   obj   - Initialized instance of the class

            obj.array = array;
        end
        
        function [thetas, phis] = estimate(obj, X, K, theta_grid, phi_grid, coherent)

            % Description: Estimates the 2D DOAs (Zenith, Azimuth) using the MUSIC algorithm.
            %
            % Inputs:
            %   X          - Received signal matrix (Sensors x Snapshots)
            %   K          - Expected number of signal sources
            %   theta_grid - Evaluation grid for elevation angles
            %   phi_grid   - Evaluation grid for azimuth angles
            %   coherent   - Boolean indicating if spatial smoothing is required
            %
            % Outputs:
            %   thetas     - Estimated elevation angles
            %   phis       - Estimated azimuth angles

            if nargin < 4 || nargin < 5
                error('Both theta_grid and phi_grid are required for 2D estimation.');
            end
            
            % Compute the 2D spectrum
            spectrum_2d = obj.musicSpectrumVectorized(X, K, theta_grid, phi_grid, coherent);
            
            % Locate regional maxima
            binary_mask = imregionalmax(spectrum_2d);
            [row, col] = find(binary_mask);
            
            % Extract and sort peak values
            peak_values = spectrum_2d(sub2ind(size(spectrum_2d), row, col));
            [peak_values_sorted, idx] = sort(peak_values, 'descend');
            row_sorted = row(idx);
            col_sorted = col(idx);
            
            % Filter peaks based on a relative power threshold
            threshold_db = -35; % dB
            peak_values_sorted_db = 10 * log10(peak_values_sorted / max(peak_values_sorted));
            
            row_filtered = [];
            col_filtered = [];
            
            for k = 1:K
                if peak_values_sorted_db(k) < threshold_db
                    fprintf('Warning: Peak at (%d,%d) is more than %d dB below the maximum and will be ignored.\n', ...
                        row_sorted(k), col_sorted(k), abs(threshold_db));
                    continue
                end
                row_filtered(end+1) = row_sorted(k);
                col_filtered(end+1) = col_sorted(k);
            end
            
            % Extract up to K valid peaks
            k_max = min(K, numel(row_filtered));
            row_final = row_filtered(1:k_max);
            col_final = col_filtered(1:k_max);
            
            % Map matrix indices to corresponding physical angles
            thetas = theta_grid(col_final);
            phis   = phi_grid(row_final);
        end
        
        function spectrum = musicSpectrum(obj, X, K, theta_grid, phi_grid, coherent)

            % Description: Computes the 2D MUSIC spectrum iteratively (loop-based).
            %
            % Inputs:
            %   X          - Signal matrix
            %   K          - Number of sources
            %   theta_grid - Vector of evaluation angles for elevation
            %   phi_grid   - Vector of evaluation angles for azimuth
            %   coherent   - Boolean indicating spatial smoothing necessity
            %
            % Outputs:
            %   spectrum   - Computed 2D spectrum matrix

            M = obj.array.M;
            N = obj.array.N;
            Lx = floor(M/2);
            Ly = floor(N/2);
            
            U_n = obj.calculateNoiseSubspace(X, K, M, N, Lx, Ly, coherent);
            U_n_H = U_n';
            
            spectrum = zeros(length(phi_grid), length(theta_grid));
            
            for j = 1:length(theta_grid) % Columns (theta)
                for i = 1:length(phi_grid) % Rows (phi)
                    if coherent
                        a_full = obj.array.steeringVector(angle_grid(i));
                        a = a_full(1:Lx, 1:Ly);
                    else
                        a = obj.array.steeringVector(theta_grid(j), phi_grid(i));
                    end
                    projection = (a' * U_n) * (U_n_H * a);
                    spectrum(i, j) = 1 / (abs(projection) + 1e-12);
                end
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
                R = MUSIC.spatialSmoothing2D(X, M_rows, N_cols, Lx, Ly);
                
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
    end
    
    methods (Static, Access = private)
        function U_n = calculateNoiseSubspace(X, K, Nx, Ny, Lx, Ly, coherent)

            % Description: Extracts the noise subspace from the signal correlation matrix.
            %
            % Inputs:
            %   X        - Signal matrix
            %   K        - Number of target sources
            %   Nx, Ny   - Total array dimensions
            %   Lx, Ly   - Sub-array dimensions for smoothing
            %   coherent - Boolean determining if spatial smoothing is applied
            %
            % Outputs:
            %   U_n      - Noise subspace eigenvector matrix

            if coherent
                R = MUSIC.spatialSmoothing2D(X, Nx, Ny, Lx, Ly);
            else
                R = (X * X') / size(X, 2);
            end
            [V, D] = eig(R);
            [~, idx] = sort(diag(D), 'descend');
            V = V(:, idx);
            
            U_n = V(:, K+1:end);
        end
        
        function R_fb = forwardBackwardAveraging(R)

            % Description: Applies Forward-Backward Averaging suitable for 
            %              Uniform Linear Arrays (ULA).
            %
            % Inputs:
            %   R    - Covariance matrix
            %
            % Outputs:
            %   R_fb - Forward-backward averaged covariance matrix

            M = size(R, 1);
            J = fliplr(eye(M));
            R_fb = 0.5 * (R + J * conj(R) * J);
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
        
            % Dimensionality checks
            if size(X,1) ~= Nant
                error('Inconsistent dimension: Matrix X must have Nx*Ny rows.');
            end
        
            % Compute grid limits for overlapping sub-arrays
            num_sub_x = Nx - Lx + 1;
            num_sub_y = Ny - Ly + 1;
            num_subarrays = num_sub_x * num_sub_y;
        
            Rss = zeros(Lx * Ly, Lx * Ly);
        
            % Iterate through each sub-array grouping
            for ix = 1:num_sub_x
                for iy = 1:num_sub_y
                    
                    % Generate sequential indices targeting the larger array
                    idx = [];
                    for x = ix:(ix + Lx - 1)
                        base = (x-1) * Ny;
                        idx = [idx, base + iy : base + iy + Ly - 1];
                    end
        
                    % Extract sub-array data and compute localized covariance
                    Xsub = X(idx, :);                  
                    Rsub = (Xsub * Xsub') / snaps;     
                    
                    Rss = Rss + Rsub;
                end
            end
        
            % Finalize mean aggregation
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
            ix = 1; % Reference X sub-array origin
            iy = 1; % Reference Y sub-array origin
            
            for x = ix:(ix+Lx-1)
                base = (x-1) * Ny;
                indices = [indices, base + iy : base + iy + Ly - 1];
            end
        end
    end
end