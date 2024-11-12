% Constants
M = 512;                       % Grid size
Mask = ones(M, M);             % Mask (modifiable)
L0 = 1000;                     % Outer scale (100 km) in m
l0 = 0.01;                     % Inner scale (1 cm)
rho_o = 0.05;                  % Effective coherence length (5 cm)
Cn = 1e-7;                     % Refractive index structure constant
delta = 0.001;                 % Grid spacing (cm)
wavelength = 5e-7;             % Wavelength (cm)
delta_x = delta;               % Grid spacing in x-direction
delta_y = delta;               % Grid spacing in y-direction
delta_z = 64;                  % Grid spacing in z-direction

% Calculate rho points
rho_points = 1:M/2;            % Only the first half of the grid

% Derived constants
kappa_l = 3.3 / l0;
kappa_0 = 2 * pi / L0;

% Function for fn(x)
fn = @(x) 1 + 1.802 * x - 0.254 * x.^(7/6);

% Phin function (element-wise operations)
phin = @(kappa) 0.033 * Cn^2 * fn(kappa / kappa_l) .* exp(-kappa.^2 / kappa_l^2) ./ (kappa.^2 + kappa_0^2).^(11/6);

% Generate phase screens
phz1 = real(phasescreen2(M, delta, L0, l0, wavelength, Cn, delta_x, delta_y, delta_z));
D1_x_y = str_fcn2_ft(phz1, Mask, delta);

phz2 = real(phasescreen(M, delta, L0, l0, wavelength, Cn, delta_z));
D2_x_y = str_fcn2_ft(phz2, Mask, delta);

% Extract diagonal elements using diag function for better performance
D1_rho = diag(D1_x_y(1:M/2, 1:M/2));  % Diagonal of upper-left M/2 x M/2 submatrix
D2_rho = diag(D2_x_y(1:M/2, 1:M/2));  % Diagonal of upper-left M/2 x M/2 submatrix

% Integration limits and grid setup
kappa_rho_max = 100;           % Upper limit for integration
kappa_rho_points = 1000;       % Number of points for integration

% Initialize D_theta_rho
D_theta_rho = zeros(size(rho_points));

% Precompute constants
constant_term = 8 * pi^2 * (2 * pi / wavelength)^2 * delta_z;

% Numerical integration using trapezoidal rule for each rho
for i = 1:length(rho_points)
    rho = rho_points(i) * delta;
    kappa_rho_values = linspace(0, kappa_rho_max, kappa_rho_points);
    
    % Integrand function
    integrand = @(kappa_rho) constant_term * kappa_rho .* phin(kappa_rho) .* (1 - besselj(0, rho * kappa_rho));
    
    % Perform the integration
    D_theta_rho(i) = integral(integrand, 0, kappa_rho_max); 
end

% Plot all three data sets on the same figure with logarithmic y-axis
figure;
loglog(rho_points * delta, D_theta_rho, 'LineWidth', 1.5, 'DisplayName', 'D_\theta(\rho)');
hold on;
loglog(rho_points * delta, D1_rho, 'LineWidth', 1.5, 'DisplayName', 'D1 trad');
loglog(rho_points * delta, D2_rho, 'LineWidth', 1.5, 'DisplayName', 'D2 rand');
hold off;

% Customize the plot
xlabel('\rho (cm)');
ylabel('D(\rho)');
title('Combined Plot of D_\theta(\rho), D1_\theta(\rho), and D2_\theta(\rho) with Logarithmic Scale');
legend('show');
grid on;
