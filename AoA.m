% Constants (adjust as necessary)
lambda = 5e-7; % Wavelength in meters (example value)
k = 2 * pi / lambda; % Wave number
L = 2000; % Propagation distance in meters
f = sqrt(lambda * L); % Fresnel length in meters
D = 0.4; % Aperture diameter in meters
l_0 = 0.01; % Inner scale of turbulence in meters
L_0 = 1000; % Outer scale of turbulence in meters

% Define the Airy function A(x)
A = @(x) (2 * besselj(1, x) ./ x).^2; % J1 is the first-order Bessel function of the first kind

% Phi_n function
function phi_n = Phi_n(k, l_0, L_0, Cn2)
    k_l = 3.3 / l_0; % Scaling constant for the inner scale
    k_0 = 2 * pi / L_0; % Scaling constant for the outer scale

    % Define the scaling function f_n(k/k_l)
    fn = @(x) 1 + 1.802 * x - 0.254 * x.^(7/6); % Empirical function

    % Calculate the power spectral density
    phi_n = 0.033 * Cn2 * fn(k / k_l) .* exp(-k.^2 / k_l^2) ./ (k.^2 + k_0^2).^(11/6);
end

% Range of Cn^2 values for Rytov variance calculations
Cn2_values = linspace(1e-25, 1e-17, 2); % Adjust based on realistic conditions
rytov_variance_values = 1.23 * Cn2_values .* k.^(7/6) .* L.^(11/6);

% Preallocate arrays for AoA variance calculations
aoa_variance_values = zeros(size(Cn2_values));

% Numerical integration for each Cn^2 value
for i = 1:length(Cn2_values)
    Cn2 = Cn2_values(i);
    
    % Define the integrand for the AoA variance calculation
    integrand = @(kappa) kappa.^3 .* Phi_n(kappa, l_0, L_0, Cn2) .* ...
                         (1 + (2 * pi) ./ (kappa * f).^2) .* ...
                         sin((kappa * f).^2 / (2 * pi)) .* A(D * kappa / 2);
    
    % Numerical integration with a finite upper limit
    aoa_variance_values(i) = pi^2 * L * integral(integrand, 0, 1000); % 1000 as a practical upper limit
end

% Plot AoA variance vs Rytov variance
figure;
loglog(rytov_variance_values, aoa_variance_values, 'b', 'LineWidth', 1.5);
hold on; % Maintain current plot for additional data

% Parameters for the phase screen method
M = 1024; % Grid size for phase screen
delta = 0.001; % Grid spacing in meters
delta_x = delta;
delta_y = delta;
delta_z = L / 13; % Step size for phase screen propagation (equal segments)

% Preallocate arrays for phase screen AoA variance calculations
aoa_variance_phasescreen = zeros(size(Cn2_values));
aoa_variance_phasescreen2 = zeros(size(Cn2_values));

% Number of phase screens
num_screens = 13;
screen_distance = L / num_screens; % Distance between phase screens

% Calculate AoA variance using the phase screen method
for i = 1:length(Cn2_values)
    Cn2 = Cn2_values(i);
    
    % Initialize the complex fields representing the wavefronts
    field1 = exp(1i * zeros(M, M)); % Initial wavefront with zero phase
    field2 = exp(1i * zeros(M, M)); % Initial wavefront with zero phase
    
    % Apply phase screens over the propagation distance
    for j = 1:num_screens
        phase_screen_segment = real(phasescreen(M, delta, L_0, l_0, lambda, Cn2,delta_z)); 
        phase_screen_segment2 = real(phasescreen2(M, delta, L_0, l_0, lambda, Cn2,delta_x, delta_y, delta_z));
        
        % Apply the phase screens to the fields
        field1 = field1 .* exp(1i * phase_screen_segment);
        field2 = field2 .* exp(1i * phase_screen_segment2);
        
        % Propagate the fields to the next screen
        field1 = propagate_field(field1, lambda, delta, delta_z);
        field2 = propagate_field(field2, lambda, delta, delta_z);
    end
    
    % Extract the phase from the final fields
    final_phase1 = angle(field1);
    final_phase2 = angle(field2);
    
    % Calculate gradients for AoA variance
    [grad_x1, grad_y1] = gradient(final_phase1, delta);
    [grad_x2, grad_y2] = gradient(final_phase2, delta);
    
    aoa_variance_phasescreen(i) = mean(grad_x1.^2 + grad_y1.^2, 'all');
    aoa_variance_phasescreen2(i) = mean(grad_x2.^2 + grad_y2.^2, 'all');
end

% Plot results of the phase screen methods
loglog(rytov_variance_values, aoa_variance_phasescreen, 'r', 'LineWidth', 1.5);
loglog(rytov_variance_values, aoa_variance_phasescreen2, 'g', 'LineWidth', 1.5);

% Final plot adjustments
xlabel('Rytov Variance \sigma_R^2');
ylabel('AoA Variance');
title('AoA Variance vs Rytov Variance');
legend('Rytov Method', 'Phase Screen Method traditional', 'Phase Screen Method randomized');
grid on;
hold off;

% Propagation function for wavefront
function field_out = propagate_field(field, lambda, delta, distance)
    [M, ~] = size(field);
    fx = (-M/2:M/2-1) / (M * delta); % Frequency space
    [FX, FY] = meshgrid(fx, fx);
    H = exp(-1i * pi * lambda * distance * (FX.^2 + FY.^2)); % Transfer function
    field_out = ifftshift(ifft2(fft2(fftshift(field)) .* H)); % Apply transfer function
end
