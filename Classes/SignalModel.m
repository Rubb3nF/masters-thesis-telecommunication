classdef SignalModel

    % Description: Narrowband signal model for DOA estimation.
    %
    % The model assumes narrowband signals impinging on a uniform array.
    %
    % Signal model equation:
    % x(t) = A(θ)s(t) + n(t)
    %
    % where:
    %   x(t) : M x 1 received signal vector
    %   A    : M x K array manifold matrix
    %   s(t) : K x 1 source signal vector
    %   n(t) : M x 1 additive noise vector
    
    properties
        array          % Antenna array configuration object
        num_elements   % Total number of sensors (M * N)
    end
    
    methods
        function obj = SignalModel(array)

            % Description: Constructor for the SignalModel class.
            %
            % Inputs:
            %   array - Initialized array object (e.g., UniformPlanarArray)
            %
            % Outputs:
            %   obj   - Initialized SignalModel instance

            obj.array = array;
            obj.num_elements = array.num_elements;
        end
        
        function [X, S, N] = generateSignals(obj, thetas, phis, N_snapshots, snr_db, signal_type, correlation, noise_type, seed)
            
            % Description: Generates synthetic array data with known DOAs.
            %
            % Inputs:
            %   thetas      - Zenith angles [1xK] (radians)
            %   phis        - Azimuth angles [1xK] (radians)
            %   N_snapshots - Number of temporal snapshots to generate
            %   snr_db      - Signal to Noise Ratio in dB
            %   signal_type - 'complex_sinusoid', 'complex_sinusoid_coherent', 
            %                 'random' (default), 'wideband_continuous', or 'chirp_continuous'
            %   correlation - Correlation coefficient between sources [0, 1]
            %   noise_type  - 'white' (default) or 'colored'
            %   seed        - (Optional) Random number generator seed for reproducibility
            %
            % Outputs:
            %   X - Received data matrix of size (M*N) x N_snapshots
            %   S - Source signal matrix of size K x N_snapshots
            %   N - Additive noise matrix of size (M*N) x N_snapshots

            if nargin < 4 || isempty(snr_db), snr_db = 10; end
            if nargin < 5 || isempty(signal_type), signal_type = 'complex_sinusoid'; end
            if nargin < 6 || isempty(correlation), correlation = 0.0; end
            if nargin < 7 || isempty(noise_type), noise_type = 'white'; end
            if nargin == 8 && ~isempty(seed)
                rng(seed);
            end
            
            thetas = thetas(:).';
            phis = phis(:).';
            K = numel(thetas);
            
            % Compute the array manifold matrix
            A = obj.array.steeringVector(thetas, phis);
            
            % Generate base source signals
            S = obj.generateSourceSignals(K, N_snapshots, signal_type, correlation);
            
            % Apply SNR scaling (expand scalar to vector if necessary)
            if isscalar(snr_db)
                snr_db = repmat(snr_db, 1, K);
            end
            
            for k = 1:K
                snr_linear = 10^(snr_db(k)/10);
                S(k, :) = S(k, :) * sqrt(snr_linear);
            end
            
            % Generate environmental noise
            N = obj.generateNoise(N_snapshots, noise_type);
            
            % Construct the final received signal matrix
            X = A * S + N;
        end
        
        function S = generateSourceSignals(obj, K, N, signal_type, correlation)

            % Description: Generates base waveforms for the signal sources.
            %
            % Inputs:
            %   K           - Number of sources
            %   N           - Number of temporal snapshots
            %   signal_type - Selected waveform type
            %   correlation - Desired correlation coefficient
            %
            % Outputs:
            %   S           - Base source signal matrix (K x N)
            
            switch signal_type
                case 'complex_sinusoid'
                    f_k = rand(1, K) * 0.5;
                    phases = 2 * pi * rand(1, K);
                    S = zeros(K, N);
                    for k = 1:K
                        S(k, :) = exp(1j * (2*pi*f_k(k)*(0:N-1) + phases(k)));
                    end
                    
                case 'complex_sinusoid_coherent'
                    phases = 2 * pi * rand(1, K);
                    S = zeros(K, N);
                    for k = 1:K
                        % Simple sinusoid with a fixed normalized frequency of 0.1
                        S(k, :) = exp(1j * (2*pi*0.1*(0:N-1) + phases(k)));
                    end
                    
                case 'random'
                    % Normalized Complex White Gaussian Noise (unit power)
                    S = (randn(K, N) + 1j*randn(K, N)) / sqrt(2);
                    
                case 'wideband_continuous' 
                    % Generate complex Gaussian noise (unit power)
                    white_noise = (randn(K, N) + 1j*randn(K, N)) / sqrt(2);
                    
                    % FIR Low-Pass Filter Design (Narrow Transition)
                    N_filter = 64;           % Filter order (higher order yields sharper transition)
                    Fc = 0.3;                % Normalized cutoff frequency (0 < Fc < 1)
                    h_filter = fir1(N_filter, Fc); 
                    
                    S = zeros(K, N);
                    for k = 1:K
                        % Apply filter via convolution
                        S(k, :) = conv(white_noise(k, :), h_filter, 'same'); 
                        
                        % Normalize filtered signal power
                        P_k = mean(abs(S(k, :)).^2);
                        if P_k > 0
                             S(k, :) = S(k, :) / sqrt(P_k);
                        end
                    end
                    
                case 'chirp_continuous'
                    T_total = N;
                    Fs = 1;                 
                    t = (0:N-1) / Fs;      
                    S = zeros(K, N);
                    
                    % Variation ranges to induce incoherence
                    B_base = 0.5 * Fs;          % Base bandwidth
                    B_variation = 0.1 * Fs;     % Bandwidth variation range 
                    F0_max = 0.2 * Fs;          % Maximum initial frequency
                    
                    for k = 1:K
                        % Random initial frequency (f_start)
                        f_start = rand() * F0_max; 
                        
                        % B_k varies around the base bandwidth
                        B_k = B_base + (rand() - 0.5) * B_variation; 
                        k_chirp = B_k / T_total; % Unique chirp rate for source k
                        
                        s_chirp_continuous = exp(1j * (2*pi*f_start*t + pi * k_chirp * t.^2));
                        
                        % Power normalization
                        P_sig = mean(abs(s_chirp_continuous).^2);
                        if P_sig > 0
                            s_chirp_continuous = s_chirp_continuous / sqrt(P_sig);
                        end
                    
                        phase_offset = exp(1j * 2*pi*rand(1));
                        S(k, :) = s_chirp_continuous * phase_offset;
                    end
                    
                otherwise
                    error('Unknown signal type: %s', signal_type);
            end
            
            % Inject correlation among sources if specified
            if correlation > 0 && K > 1
                rho = correlation;
                for k = 2:K
                    independent = (randn(1, N) + 1j*randn(1, N)) / sqrt(2);
                    S(k, :) = rho * S(1, :) + sqrt(1 - rho^2) * independent;
                    
                    % Re-normalize power after mixing
                    P_k = mean(abs(S(k, :)).^2);
                    if P_k > 0
                         S(k, :) = S(k, :) / sqrt(P_k);
                    end
                end
            end
        end
        
        function N = generateNoise(obj, N_snapshots, noise_type)

            % Description: Generates additive sensor noise.
            %
            % Inputs:
            %   N_snapshots - Number of temporal snapshots
            %   noise_type  - Target noise profile ('white' or 'colored')
            %
            % Outputs:
            %   N           - Additive noise matrix (M*N x N_snapshots)

            switch noise_type
                case 'white'
                    N = (randn(obj.num_elements, N_snapshots) + 1j*randn(obj.num_elements, N_snapshots)) / sqrt(2);
                case 'colored'
                    white_noise = (randn(obj.num_elements, N_snapshots) + 1j*randn(obj.num_elements, N_snapshots)) / sqrt(2);
                    h = [1, 0.5];
                    N = zeros(size(white_noise));
                    for m = 1:obj.num_elements
                        N(m, :) = conv(white_noise(m, :), h, 'same');
                    end
                otherwise
                    error('Unknown noise type: %s', noise_type);
            end
        end
        
        function R = sampleCovariance(obj, X)

            % Description: Computes the sample covariance matrix from received data.
            %
            % Inputs:
            %   X - Received data matrix (M*N x Snapshots)
            %
            % Outputs:
            %   R - Sample covariance matrix (M*N x M*N)

            R = (X * X') / size(X, 2);
        end
        
        function R = theoreticalCovariance(obj, thetas, phis, powers, noise_power)

            % Description: Computes the theoretical (asymptotic) covariance matrix.
            %
            % Inputs:
            %   thetas      - Zenith angles [1xK] (radians)
            %   phis        - Azimuth angles [1xK] (radians)
            %   powers      - Vector of source signal powers
            %   noise_power - Scalar variance of the noise (default: 1.0)
            %
            % Outputs:
            %   R           - Theoretical covariance matrix (M*N x M*N)

            if nargin < 5, noise_power = 1.0; end
            
            A = obj.array.steeringVector(thetas, phis);
            P = diag(powers);
            
            R = A * P * A' + noise_power * eye(obj.num_elements);
        end
        
        function snr_db = snrEstimate(obj, X, thetas, phis)

            % Description: Estimates the empirical Signal-to-Noise Ratio.
            %
            % Inputs:
            %   X      - Received data matrix (M*N x Snapshots)
            %   thetas - Known or estimated zenith DOAs (radians)
            %   phis   - Known or estimated azimuth DOAs (radians)
            %
            % Outputs:
            %   snr_db - Vector containing estimated SNR in dB for each source

            R = obj.sampleCovariance(X);
            eigenvals = sort(real(eig(R)), 'descend');
            K = numel(thetas); 
            
            if K + 1 > length(eigenvals)
                % Insufficient degrees of freedom to estimate noise subspace
                noise_power = 1e-10; 
            else
                noise_power = mean(eigenvals(K+1:end));
            end
            
            signal_power = mean(eigenvals(1:K));
            
            snr_linear = (signal_power - noise_power) / noise_power;
            snr_linear = max(snr_linear, 1e-10); % Bound to prevent log(0)
            
            snr_db = repmat(10 * log10(snr_linear), 1, K);
        end
    end
end