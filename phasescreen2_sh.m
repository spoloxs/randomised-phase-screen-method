function theta_H = phasescreen2_sh(M, delta, L0, l0, wavelength, Cn, delta_x, delta_y, delta_z, N_p)
    % Compute theta_R(j, l) using the existing phasescreen2 function
    theta_R = phasescreen2(M, delta, L0, l0, wavelength, Cn, delta_x, delta_y, delta_z);
    
    % Set up additional parameters for theta_out and theta_in
    del_k = 2 * pi / (M * delta);            % Frequency grid spacing [1/m]
    kx = (-M/2 : M/2-1) * del_k;             % Frequency grid [1/m]
    [kx, ky] = meshgrid(kx, kx);             % Create a 2D frequency grid
    
    % Compute the theta_out(j, l, p) and theta_in(j, l)
    theta_out = zeros(M, M);
    theta_in = zeros(M, M);
    
    % Compute theta_out
    for p = 1:N_p
        scale = 3^(-p);
        for n = -1:1
            for m = -1:1
                if ~(n == 0 && m == 0) % Exclude the delta function at (0, 0)
                    delta_kx = (randn - 0.5) * del_k; % Random phase offsets
                    delta_ky = delta_kx; % Same random offset for ky
                    theta_out = theta_out + scale * ...
                        computeCtilde(kx, ky, M, l0, L0, Cn, delta_kx, delta_ky, p, delta_x, delta_y);
                end
            end
        end
    end
    
    % Compute theta_in
    scale_in = 3^(-(N_p + 1));
    for n = -1:1
        for m = -1:1
            if ~(n == 0 && m == 0) % Exclude the delta function at (0, 0)
                delta_kx = (randn - 0.5) * del_k; % Random phase offsets
                delta_ky = delta_kx; % Same random offset for ky
                theta_in = theta_in + scale_in * ...
                    computeCtilde(kx, ky, M, l0, L0, Cn, delta_kx, delta_ky, N_p + 1, delta_x, delta_y);
            end
        end
    end
    
    % Final summation to get theta_H
    theta_H = theta_R + theta_out + theta_in;
end

function C_tilde_val = computeCtilde(kx, ky, M, l0, L0, Cn, delta_kx, delta_ky, p, delta_x, delta_y)
    % Helper function to compute the C_tilde term for theta_out and theta_in
    a = randn(M, M); % Real part (Gaussian distributed)
    b = randn(M, M); % Imaginary part (Gaussian distributed)
    z = a + 1i * b;
    C_tilde = z .* sqrt(2 * pi * (2 * pi / (M * delta_x))^2 * computePhi_n(kx + delta_kx / 3^p, ky + delta_ky / 3^p, l0, L0, Cn));
    C_tilde_val = M^2 * ifft2(C_tilde) .* exp(1i * ((0:M-1).' * delta_x .* (kx + delta_kx / 3^p) + (0:M-1).' * delta_y .* (ky + delta_ky / 3^p)));
end
