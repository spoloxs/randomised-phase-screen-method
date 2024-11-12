function phz = phasescreen2(M, delta, L0, l0, wavelength, Cn, delta_x, delta_y, delta_z)
    % Set up the frequency grid
    del_k = 2 * pi / (M * delta);            % Frequency grid spacing [1/m]
    kx = (-M/2 : M/2-1) * del_k;             % Frequency grid [1/m]
    [kx, ky] = meshgrid(kx, kx);             % Create a 2D frequency grid

    constant_value_kx = (rand) * del_k/4;  % Generate a single random value
    delta_kx = constant_value_kx * ones(M, M);  % Assign the constant value to all points
    delta_ky = delta_kx;  % Assign theant value to all points
    
    % Generate Gaussian distributed random values for the real and imaginary parts
    a = randn(M, M);  % Real part
    b = randn(M, M);  % Imaginary part
    
    % Combine 'a' and 'b' to form the complex matrix 'z'
    z = a + 1i * b;
    
    % Precompute the common part of the equation for C_tilde
    kx_delta_kx = kx + delta_kx;  % Add delta_kx to each kx value
    ky_delta_ky = ky + delta_ky;  % Add delta_ky to each ky value
    
    % Precompute the value of computePhi_n for each (kx, ky)
    Phi_n_vals = computePhi_n(kx_delta_kx, ky_delta_ky, l0, L0, Cn);
    
    % Calculate C_tilde in one step by vectorizing the operations
    C_tilde = z .* sqrt(2 * pi * delta_z * del_k^2 * Phi_n_vals) * (2 * pi / wavelength);
    
    % Compute the inverse FFT to get the phase screen in the spatial domain
    C = M.^2 * ifft2(C_tilde);
    
    % Precompute the indices for theta_R
    [m, n] = size(C);  % Get the size of matrix C
    [X, Y] = meshgrid(0:m-1, 0:n-1);  % Create mesh grids for the indices
    theta_R = exp(1i * (X * delta_x .* delta_kx + Y * delta_y .* delta_ky)) .* C;

    % Return the phase screen (optional: theta_R or C can be returned)
    phz = theta_R;
end
