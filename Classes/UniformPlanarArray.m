classdef UniformPlanarArray

    % Description: Uniform Planar Array (UPA) class for Direction of Arrival 
    %              estimation.
    %
    % The general steering vector formula is implemented as:
    % a(theta, phi) = exp(j * k' * P)
    %
    % Angle Conventions:
    %   theta : Zenith angle (from Z-axis, 0 to pi)
    %   phi   : Azimuth angle (in X-Y plane, from X-axis, 0 to 2pi)

    
    properties
        M                 % Number of elements in the x-direction
        N                 % Number of elements in the y-direction
        dx                % Inter-element spacing in x (wavelengths)
        dy                % Inter-element spacing in y (wavelengths)
        fc                % Carrier frequency (Hz)
        c                 % Speed of light (m/s)
        wavelength        % Signal wavelength (m)
        num_elements      % Total number of antenna elements (M * N)
        
        % Element positions matrix P
        % Dimension: 3 x (M*N) matrix formatted as [x_pos; y_pos; z_pos]
        % Units: Meters
        element_positions 
    end
    
    methods
        function obj = UniformPlanarArray(M, N, dx, dy, fc, c)

            % Description: Constructor for the UniformPlanarArray class.
            %
            % Inputs:
            %   M  - Number of elements in the x-direction
            %   N  - Number of elements in the y-direction
            %   dx - Inter-element spacing in x (wavelengths, default: 0.5)
            %   dy - Inter-element spacing in y (wavelengths, default: dx)
            %   fc - Carrier frequency in Hz (default: 1e9)
            %   c  - Speed of light in m/s (default: 3e8)
            %
            % Outputs:
            %   obj - Initialized UniformPlanarArray instance

            if nargin < 2
                error('Dimensions M and N must be specified.'); 
            end
            if nargin < 3, dx = 0.5; end
            if nargin < 4, dy = dx;  end 
            if nargin < 5, fc = 1e9; end
            if nargin < 6, c = 3e8;  end
            
            obj.M = M;
            obj.N = N;
            obj.dx = dx;
            obj.dy = dy;
            obj.fc = fc;
            obj.c = c;
            obj.wavelength = c / fc;
            obj.num_elements = M * N;
            
            % Generate spatial grid for elements
            [m_grid, n_grid] = ndgrid(0:(M-1), 0:(N-1));
            
            % Compute physical coordinates (1 x MN vectors)
            x_pos = (m_grid(:).' * obj.dx) * obj.wavelength;
            y_pos = (n_grid(:).' * obj.dy) * obj.wavelength;
            z_pos = zeros(1, obj.num_elements);
            
            % Store positions in meters
            obj.element_positions = [x_pos; y_pos; z_pos];
        end
        
        function a = steeringVector(obj, theta, phi)

            % Description: Computes the array steering vector(s) for given directions.
            %
            % The wave vector k is defined as:
            % $$\mathbf{k} = \frac{2\pi}{\lambda} [\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta]^T$$
            %
            % Inputs:
            %   theta - Zenith angle(s) from z-axis [1xK] (radians)
            %   phi   - Azimuth angle(s) in x-y plane [1xK] (radians)
            %
            % Outputs:
            %   a     - Computed steering matrix of size (M*N) x K

            theta = theta(:).'; % Enforce row vector formulation (1xK)
            phi = phi(:).';     % Enforce row vector formulation (1xK)
            
            if length(theta) ~= length(phi)
                error('Theta and phi arrays must contain the same number of elements.');
            end
            
            % Wave vector components for all K directions
            k_magnitude = 2 * pi / obj.wavelength;
            
            k_x = k_magnitude * sin(theta) .* cos(phi);
            k_y = k_magnitude * sin(theta) .* sin(phi);
            k_z = k_magnitude * cos(theta);
            
            % Aggregation into a (3 x K) matrix
            k_vectors = [k_x; k_y; k_z];
            
            % Calculation of the spatial phase matrix via dot product
            % Dimension mapping: (K x 3) * (3 x MN) -> (K x MN)
            phase = k_vectors.' * obj.element_positions;
            
            % Evaluation of the complex exponential
            a = exp(1j * phase.');
        end
        
        function A = arrayManifold(obj, theta_list, phi_list)

            % Description: Convenience wrapper to compute the array manifold matrix 
            %              over a set of evaluated angles.
            %
            % Inputs:
            %   theta_list - Zenith angle(s) [1xK] (radians)
            %   phi_list   - Azimuth angle(s) [1xK] (radians)
            %
            % Outputs:
            %   A          - Array manifold matrix of size (M*N) x K

            A = obj.steeringVector(theta_list, phi_list);
        end
        
        function [u, v] = spatialSines(obj, theta, phi)

            % Description: Converts standard physical DOA to spatial sines.
            %
            % Mathematical mapping:
            % $$u = \sin(\theta)\cos(\phi)$$
            % $$v = \sin(\theta)\sin(\phi)$$
            %
            % Inputs:
            %   theta - Zenith angle (radians)
            %   phi   - Azimuth angle (radians)
            %
            % Outputs:
            %   u, v  - Spatial sine domain coordinates

            u = sin(theta) .* cos(phi);
            v = sin(theta) .* sin(phi);
        end
        
        function [theta, phi] = anglesFromSines(obj, u, v)

            % Description: Transforms spatial sines back to physical DOA angles.
            %
            % Inputs:
            %   u, v  - Spatial sine domain coordinates
            %
            % Outputs:
            %   theta - Reconstructed zenith angle (radians)
            %   phi   - Reconstructed azimuth angle (radians)

            
            % Zenith angle
            sin_theta_sq = u.^2 + v.^2;
            
            % Hard clipping to mitigate numerical/noise artifacts (>1)
            sin_theta_sq(sin_theta_sq > 1) = 1;
            
            theta = asin(sqrt(sin_theta_sq));
            
            % Azimuth angle
            phi = atan2(v, u);
        end
        
        function info = get_info(obj)

            % Description: Aggregates and returns core array configuration parameters.
            %
            % Outputs:
            %   info - Struct containing dimensional and physical array properties

            info = struct( ...
                'M', obj.M, ...
                'N', obj.N, ...
                'num_elements', obj.num_elements, ...
                'dx_wavelengths', obj.dx, ...
                'dy_wavelengths', obj.dy, ...
                'dx_meters', obj.dx * obj.wavelength, ...
                'dy_meters', obj.dy * obj.wavelength, ...
                'fc', obj.fc, ...
                'wavelength', obj.wavelength, ...
                'element_positions_meters', obj.element_positions ...
            );
        end
    end
end