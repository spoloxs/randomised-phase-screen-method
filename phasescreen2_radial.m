function phz = phasescreen2_radial(M, delta, L0, l0, wavelength, Cn, delta_x, delta_y, delta_z)
    % Set up the radial frequency grid
    del_k = 2 * pi / (M * delta);          % Frequency grid spacing [1/m]
    kx = (-M/2 : M/2-1) * del_k;           % Frequency grid [1/m]
    [kx, ky] = meshgrid(kx, kx);           % Create a 2D frequency grid

    % Compute the radial frequency grid
    kr = sqrt(kx.^2 + ky.^2);              % Radial frequency components (magnitude)
    theta_r = atan2(ky, kx);               % Angular components (angle)

    % Calculate the phase screen in the frequency domain (radial coordinates)
    C_tilde = zeros(M, M);                 % Initialize phase screen in frequency domain

    % Random phase offset parameters (radial)
    constant_value_kx = (randn - 0.5) * del_k;  % Generate a single random value for the offset
    
    delta_kx = constant_value_kx * ones(M, M);  % Apply the constant value to all points
    delta_ky = delta_kx;                     % Assign the constant value to all points

    for i = 1:M
        for j = 1:M
            z = (randn + 1i * randn);   % Complex random phase
            % Compute the phase screen using radial frequencies
            C_tilde(i, j) = z * ...
                sqrt(2 * pi * delta_z * del_k^2 * computePhi_n(kr(i,j) + delta_kx(i,j), kr(i,j) + delta_ky(i,j), l0, L0, Cn)) * ...
                (2 * pi / wavelength);
        end
    end

    % Compute the inverse FFT to get the phase screen in the spatial domain
    C = M^2 * ifft2(C_tilde);             % Inverse Fourier transform to spatial domain

    % Initialize theta_R matrix (same size as C)
    [m, n] = size(C);  % Get the size of matrix C
    theta_R = zeros(m, n);  % Initialize theta_R with the same size as C
    
    % Calculate theta_R in the radial system
    for i = 1:m
        for j = 1:n
            % Calculate theta_R using radial frequency components
            theta_R(i, j) = exp(1i * (kr(i,j) * delta_x + kr(i,j) * delta_y)) * C(i, j);
        end
    end

    % Return the phase screen (optional: theta_R or C can be returned)
    phz = theta_R;

end
