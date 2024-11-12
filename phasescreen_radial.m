function phz = phasescreen_radial(M, delta, L0, l0, wavelength, Cn, delta_z)
    % Set up the radial frequency grid
    del_k = 2 * pi / (M * delta);          % Frequency grid spacing [1/m]
    kx = (-M/2 : M/2-1) * del_k;           % Frequency grid [1/m]
    [kx, ky] = meshgrid(kx, kx);           % Create a 2D frequency grid

    % Compute the radial frequency grid
    kr = sqrt(kx.^2 + ky.^2);              % Radial frequency components (magnitude)
    theta_r = atan2(ky, kx);               % Angular components (angle)

    % Calculate the phase screen in the frequency domain (radial coordinates)
    C_tilde = zeros(M, M);                 % Initialize phase screen in frequency domain

    for i = 1:M
        for j = 1:M
            z = (randn + 1i * randn);     % Complex random phase
            % Compute the phase screen using radial frequencies
            C_tilde(i, j) = z * ...
                sqrt(2 * pi * delta_z * del_k^2 * computePhi_n(kr(i,j), kr(i,j), l0, L0, Cn)) * ...
                (2 * pi / wavelength);
        end
    end

    % Compute the inverse FFT to get the phase screen in the spatial domain
    C = M^2 * ifft2(C_tilde);             % Inverse Fourier transform to spatial domain

    % Return the phase screen (optional: theta_R or C can be returned)
    phz = C;
end
