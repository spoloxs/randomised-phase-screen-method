% Define parameters
r0 = 0.15; % Fried's parameter in meters
N = 512; % Grid size
delta = 0.01; % Grid spacing in meters
L0 = 100; % Outer scale of turbulence in meters
l0 = 0.01; % Inner scale of turbulence in meters

% Generate the phase screens
[phz_lo, phz_hi] = ft_sh_phase_screen(r0, N, delta, L0, l0);

% Create a mask for the structure function calculation (e.g., circular)
mask = ones(N); % Simple full-aperture mask for demonstration

% Compute the structure function for phz_hi and phz_lo
D_hi = str_fcn2_ft(phz_hi, mask, delta);
D_lo = str_fcn2_ft(phz_lo, mask, delta);

% Create a spatial grid for plotting
x = (-N/2 : N/2-1) * delta;

% Plot D versus x for the central row using log-log scale
figure;
loglog(x, D_hi(round(N/2), :), 'b', 'DisplayName', 'High-Frequency Screen');
hold on;
xlabel('x (meters)');
ylabel('D(x) (Structure Function)');
title('Structure Function D(x) vs. x');
legend('Location', 'best');
grid on;
hold off;


function [phz_lo, phz_hi] = ft_sh_phase_screen(r0, N, delta, L0, l0)
    % Calculate screen size in meters
    D = N * delta;
    
    % High-frequency phase screen using FFT method
    phz_hi = ft_phase_screen(r0, N, delta, L0, l0);
    
    % Spatial grid [m]
    [x, y] = meshgrid((-N/2 : N/2-1) * delta);
    
    % Initialize low-frequency phase screen
    phz_lo = zeros(size(phz_hi));
    
    % Loop over frequency grids with spacing 1/(3^p * D)
    for p = 1:3
        % Frequency grid spacing [1/m]
        del_f = 1 / (3^p * D);
        fx = (-1 : 1) * del_f;
        
        % Frequency grid [1/m]
        [fx, fy] = meshgrid(fx);
        [~, f] = cart2pol(fx, fy); % Polar coordinates
        
        fm = 5.92 / l0 / (2 * pi); % Inner scale frequency [1/m]
        f0 = 1 / L0; % Outer scale frequency [1/m]
        
        % Modified von Karman atmospheric phase PSD
        PSD_phi = 0.023 * r0^(-5/3) * exp(-(f / fm).^2) ./ (f.^2 + f0^2).^(11/6);
        PSD_phi(2, 2) = 0; % Avoid singularity at the origin
        
        % Random Fourier coefficients
        cn = (randn(3) + 1i * randn(3)) .* sqrt(PSD_phi) * del_f;
        
        % Initialize subharmonics
        SH = zeros(N);
        
        % Loop over all frequencies on this grid
        for ii = 1:9
            SH = SH + cn(ii) * exp(1i * 2 * pi * (fx(ii) * x + fy(ii) * y));
        end
        
        % Accumulate subharmonics
        phz_lo = phz_lo + SH;
    end
    
    % Subtract mean to remove piston term
    phz_lo = real(phz_lo) - mean(real(phz_lo(:)));
end

function phz = ft_phase_screen(r0, N, delta, L0, l0)
    % Setup the PSD
    del_f = 1 / (N * delta); % Frequency grid spacing [1/m]
    fx = (-N/2 : N/2-1) * del_f; % Frequency grid [1/m]
    [fx, fy] = meshgrid(fx);
    [th, f] = cart2pol(fx, fy); % Polar grid

    fm = 5.92 / l0 / (2 * pi); % Inner scale frequency [1/m]
    f0 = 1 / L0; % Outer scale frequency [1/m]

    % Modified von Karman atmospheric phase PSD
    PSD_phi = 0.023 * r0^(-5/3) * exp(-(f / fm).^2) ./ (f.^2 + f0^2).^(11/6);
    PSD_phi(N/2 + 1, N/2 + 1) = 0; % Center frequency set to 0

    % Random Fourier coefficients
    cn = (randn(N) + 1i * randn(N)) .* sqrt(PSD_phi) * del_f;

    % Synthesize the phase screen
    phz = real(ifft2(cn));
end
