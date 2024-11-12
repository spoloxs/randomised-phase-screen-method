function phz = phasescreen(M, delta, L0, l0, wavelength, Cn, delta_z)
    % Set up the frequency grid
    del_k = 2 * pi / (M * delta);            % Frequency grid spacing [1/m]
    kx = (-M/2 : M/2-1) * del_k;             % Frequency grid [1/m]
    [kx, ky] = meshgrid(kx, kx);             % Create a 2D frequency grid

    % Precompute constants for the calculation
    constant_term = sqrt(2 * pi * delta_z * del_k^2) * (2 * pi / wavelength);

    % Generate Gaussian random values for the real and imaginary parts
    z = randn(M, M) + 1i * randn(M, M);  % Complex random matrix
    
    % Calculate Phi_n for the entire frequency grid (instead of looping)
    Phi_n_vals = computePhi_n(kx, ky, l0, L0, Cn);

    % Compute C_tilde in one step (using element-wise multiplication)
    C_tilde = z .* (constant_term * sqrt(Phi_n_vals));

    % Compute the inverse FFT to get the phase screen in the spatial domain
    C = M.^2 * ifft2(C_tilde);

    % Return the phase screen
    phz = C;
end
