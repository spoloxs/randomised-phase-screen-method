% Parameters
wavelength = 500e-9;       % Example wavelength [m]
delta = 1;              % Grid spacing [m]
M = 2048;                  % Grid size (NxN)
l0 = 0.01;                 % Inner scale of turbulence [m]
Cn = 1e-7;               % Structure constant [m^(-2/3)]
delta_x = delta;          % Spacing in x-direction [m]
delta_y = delta;          % Spacing in y-direction [m]
delta_z = delta;              % Propagation distance [m]

% Define the range of L0 values to test
L0_values = linspace(1000, 200000, 200);  % Example L0 range [m]
num_trials = 10;
% Initialize array to store RMS errors
errors = zeros(size(L0_values));
errors2 = zeros(size(L0_values));

% Loop through L0 values and compute RMS errors
for i = 1:length(L0_values)
    L0 = L0_values(i);  % Update L0

    % Initialize error accumulator for 10 screens
    total_ex = 0;
    total_ey = 0;
    total_ex2 = 0;
    total_ey2 = 0;
    
    % Generate and process 10 phase screens for averaging
    for screen_num = 1:num_trials
        % Generate the phase screen (assuming phz function is defined)
        phz = phasescreen(M, delta, L0, l0, wavelength, Cn, delta_z);
        phz2 = phasescreen2(M, delta, L0, l0, wavelength, Cn, delta_z, delta_y, delta_z);
        
        % Precompute D_theta(rho)
        k = 2 * pi / wavelength;  % Wavenumber
        rho = (1:M/2) * delta * sqrt(2);  % Radial distance

        % Precompute D_theta for all rho values
        D_theta = zeros(1, M/2);
        for idx = 1:length(rho)
            D_theta(idx) = computeDtheta(rho(idx), k, delta_z, l0, L0, Cn);
        end

        % Calculate the structure function for the entire phase screen
        D = str_fcn2_ft(real(phz), ones(M, M), delta);
        D2 = str_fcn2_ft(real(phz2), ones(M, M), delta);
        
        % Extract Dx and Dy from the structure function
        Dx = diag(D(1:M/2, 1:M/2));  % Diagonal of the upper-left (M/2) x (M/2) submatrix of D
        Dy = diag(D(1:M/2, 1:M/2)); % Diagonal of the upper-left (M/2) x (M/2) submatrix of D2
        Dx2 = diag(D2(1:M/2, 1:M/2));  % Adjust indexing to match size of D_theta
        Dy2 = diag(D2(1:M/2, 1:M/2)); 
        
        % Ensure Dx and Dy are of the same length as D_theta for error computation
        Dx = Dx(:);  % Reshape Dx to a column vector
        Dy = Dy(:);  % Reshape Dy to a column vector
        Dx2 = Dx2(:);  % Reshape Dx2 to a column vector
        Dy2 = Dy2(:);  % Reshape Dy2 to a column vector

        % RMS Error Metrics for current screen
        ex = 0;
        ey = 0;
        for j = 1:length(Dx)
            ex = ex + ((Dx(j) - D_theta(j)) / D_theta(j)).^2;
            ey = ey + ((Dy(j) - D_theta(j)) / D_theta(j)).^2;
        end
        ex = sqrt(ex / length(Dx));  % Correct RMS error calculation (mean instead of multiplied by 2/M)
        ey = sqrt(ey / length(Dy));  % Correct RMS error calculation
        
        % Add to the total RMS errors for averaging
        total_ex = total_ex + ex;
        total_ey = total_ey + ey;

        % RMS Error Metrics for current screen (phz2)
        ex2 = 0;
        ey2 = 0;
        for j = 1:length(Dx2)
            ex2 = ex2 + ((Dx2(j) - D_theta(j)) / D_theta(j)).^2;
            ey2 = ey2 + ((Dy2(j) - D_theta(j)) / D_theta(j)).^2;
        end
        ex2 = sqrt(ex2 / length(Dx2));  % Correct RMS error calculation
        ey2 = sqrt(ey2 / length(Dy2));  % Correct RMS error calculation
        
        % Add to the total RMS errors for averaging
        total_ex2 = total_ex2 + ex2;
        total_ey2 = total_ey2 + ey2;
    end
    
    % Average RMS error for the 10 screens
    Ex_avg = total_ex / num_trials;
    Ey_avg = total_ey / num_trials;
    
    % Combined RMS Error for the 10 screens
    E_avg = 100 * (Ex_avg + Ey_avg) / 2;

    % Store the average RMS error for this L0 value
    errors(i) = E_avg;

    % Average RMS error for the 10 screens (phz2)
    Ex_avg2 = total_ex2 / num_trials;
    Ey_avg2 = total_ey2 / num_trials;
    
    % Combined RMS Error for the 10 screens (phz2)
    E_avg2 = 100 * (Ex_avg2 + Ey_avg2) / 2;

    % Store the average RMS error for this L0 value (phz2)
    errors2(i) = E_avg2;
end

% Plot RMS Error vs L0 / (M * delta_x)
figure;
loglog(L0_values / (M * delta_x), errors, 'LineWidth', 2);
hold on;
loglog(L0_values / (M * delta_x), errors2, 'LineWidth', 2);
hold off;
xlabel('L_0 / (M * \delta_x)');
ylabel('Average RMS Error (%)');
title('Average RMS Error vs. L_0 / (M * \delta_x)');
legend('phz', 'phz2'); % Add legend for clarity
grid on;

% Helper function to compute D_theta(rho)
function D_theta = computeDtheta(rho, k, delta_z, l0, L0, Cn)
    kappa_l = 3.3 / l0;
    kappa_0 = 2 * pi / L0;
    
    integrand = @(kappa_rho) kappa_rho .* (0.033 * Cn^2 .* ...
        (1 + 1.802 * (kappa_rho / kappa_l) - 0.254 * (kappa_rho / kappa_l).^(7/6)) .* ...
        exp(-kappa_rho.^2 / kappa_l^2) ./ (kappa_rho.^2 + kappa_0^2).^(11/6)) .* ...
        (1 - besselj(0, rho * kappa_rho));
    
    D_theta = 8 * pi^2 * k^2 * delta_z * integral(integrand, 0, Inf);
end
