classdef JitteredUPA
    % Description: Jittered Uniform Planar Array class for Direction 
    %              of Arrival estimation.
    %
    % This array utilizes the complete (M x N) grid but applies a random 
    % spatial displacement (jitter) to each antenna element. This breaks 
    % the array's periodicity, helping to mitigate grating lobes when the 
    % inter-element spacing exceeds half a wavelength (d > lambda/2).
    %
    % The general steering vector formula is implemented as:
    % a(theta, phi) = exp(j * k^T * P)
    %
    % Angle Conventions:
    %   theta : Zenith angle (from Z-axis, 0 to pi)
    %   phi   : Azimuth angle (in X-Y plane, from X-axis, 0 to 2pi)
    
    properties
        M                 % Number of elements in the x-direction
        N                 % Number of elements in the y-direction
        dx                % Base inter-element spacing in x (wavelengths)
        dy                % Base inter-element spacing in y (wavelengths)
        jitter_percent    % Percentage of maximum spatial displacement (0 to 1)
        fc                % Carrier frequency (Hz)
        c                 % Speed of light (m/s)
        wavelength        % Signal wavelength (m)
        
        num_elements      % Total number of antenna elements (M * N)
        
        % Element positions matrix P (includes jitter)
        % Dimension: 3 x (M*N) matrix formatted as [x_pos; y_pos; z_pos]
        % Units: Meters
        element_positions 
    end
    
    methods
        function obj = JitteredUPA(M, N, dx, dy, jitter_percent, fc, c)

            % Description: Constructor for the JitteredUPA class.
            %
            % Inputs:
            %   M              - Number of elements in the x-direction
            %   N              - Number of elements in the y-direction
            %   dx             - Base spacing in x (wavelengths, default: 0.5)
            %   dy             - Base spacing in y (wavelengths, default: dx)
            %   jitter_percent - Allowed displacement percentage (default: 0)
            %   fc             - Carrier frequency in Hz (default: 1e9)
            %   c              - Speed of light in m/s (default: 3e8)
            %
            % Outputs:
            %   obj            - Initialized JitteredUPA instance
            
            if nargin < 2
                error('Dimensions M and N must be specified.'); 
            end
            if nargin < 3, dx = 0.5; end
            if nargin < 4, dy = dx;  end
            if nargin < 5, jitter_percent = 0; end
            if nargin < 6, fc = 1e9; end
            if nargin < 7, c = 3e8;  end
            
            obj.M = M; 
            obj.N = N;
            obj.dx = dx; 
            obj.dy = dy;
            obj.jitter_percent = jitter_percent;
            obj.fc = fc; 
            obj.c = c;
            obj.wavelength = c / fc;
            obj.num_elements = M * N;
            
            % Generate Complete Base Grid
            [m_grid, n_grid] = ndgrid(0:(M-1), 0:(N-1));
            
            % Column vectors mapping all theoretical positions
            base_x = (m_grid(:).' * obj.dx) * obj.wavelength;
            base_y = (n_grid(:).' * obj.dy) * obj.wavelength;
            
            % Apply Random Jitter to All Elements
            % Calculate maximum absolute displacement limits in meters
            max_j_x = obj.dx * obj.wavelength * obj.jitter_percent;
            max_j_y = obj.dy * obj.wavelength * obj.jitter_percent;
            
            % Generate uniform noise centered at 0: range [-max_j, max_j]
            noise_x = (rand(1, obj.num_elements) - 0.5) * 2 * max_j_x;
            noise_y = (rand(1, obj.num_elements) - 0.5) * 2 * max_j_y;
            
            % Displace elements
            final_x = base_x + noise_x;
            final_y = base_y + noise_y;
            
            % Store physical positions in meters (z-axis remains 0 for planar)
            obj.element_positions = [final_x; final_y; zeros(1, obj.num_elements)];
        end
        
        function a = steeringVector(obj, theta, phi)

            % Description: Computes the array steering vector(s) for given directions,
            %              accounting for the exact jittered positions.
            %
            % Inputs:
            %   theta - Zenith angle(s) from z-axis [1xK] (radians)
            %   phi   - Azimuth angle(s) in x-y plane [1xK] (radians)
            %
            % Outputs:
            %   a     - Computed steering matrix of size (M*N) x K

            theta = theta(:).'; 
            phi = phi(:).';
            
            if length(theta) ~= length(phi)
                error('Theta and phi arrays must contain the same number of elements.');
            end
            
            k_mag = 2 * pi / obj.wavelength;
            
            k_x = k_mag * sin(theta) .* cos(phi);
            k_y = k_mag * sin(theta) .* sin(phi);
            k_z = k_mag * cos(theta);
            
            % Exact phase calculation based on physical element positions
            % Dimension mapping: (K x 3) * (3 x MN) -> (K x MN)
            phase = [k_x; k_y; k_z].' * obj.element_positions;
            
            % Transpose required to output an (MN x K) matrix
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
        
        function plotGeometry(obj)
            figure('Color', 'w');
            plot(obj.element_positions(1,:), obj.element_positions(2,:), ...
                 '.', 'Color', 'b', 'MarkerSize', 10);
             
            title(sprintf('Full Jittered UPA (%d x %d) - Jitter %.0f%%', ...
                  obj.M, obj.N, obj.jitter_percent * 100));
            xlabel('x Position (m)'); 
            ylabel('y Position (m)'); 
            axis equal; 
            grid on;
        end
    end
end