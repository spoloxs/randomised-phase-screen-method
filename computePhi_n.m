function phi_n = computePhi_n(kx, ky, l_0, L_0, Cn)
    % Compute the scaling constants
    k_l = 3.3 / l_0;
    k_0 = 2 * pi / L_0;
    
    % Compute the modulus of the wave vector
    k_mod = sqrt(kx.^2 + ky.^2);
    % Define the scaling function f_n(k/k_l)
    fn = @(x) 1 + 1.802 * x - 0.254 * x.^(7/6);  % f_n function

    % Calculate the power spectral density
    phi_n = 0.033 * Cn^2 * fn(k_mod / k_l) .* exp(-k_mod.^2 / k_l^2) ./ (k_mod.^2 + k_0^2).^(11/6);

end
