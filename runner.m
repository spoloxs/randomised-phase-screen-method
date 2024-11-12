% Constants
N = 1024;                       % Grid size
L0 = 1000;                   % Outer scale (100 km) in cm
l0 = 0.01;                        % Inner scale (1 cm)
rho_o = 0.05;                     % Effective coherence length (5 cm)
Cn = 1e-7;                     % Refractive index structure constant
delta = 1;                 % Grid spacing (cm)
wavelength = 5e-7;             % Wavelength (cm)
delta_x = delta;               % Grid spacing in x-direction
delta_y = delta;               % Grid spacing in y-direction
delta_z = delta;                   % Grid spacing in z-direction

% % Generate the phase screen
% phz = phasescreen(N, delta, L0, l0, wavelength,Cn, delta_z);
% phz2 = phasescreen2(N, delta, L0, l0, wavelength,Cn, delta_x, delta_y, delta_z);
% % phz3 = phasescreen2_sh(N, delta, L0, l0, wavelength, Cn, delta_x, delta_y, delta_z, 1);
% % Repeat the phase screen to visualize without periodicity
% % repeated_phase_screen = repmat(real(phz), 4,4);
% % 
% % % Plot the repeated phase screen
% % figure;
% % imagesc(repeated_phase_screen);
% % colormap("jet");           % Apply the 'jet' colormap
% % colorbar; 
% % title('Repeated Randomized Phase Screen trad (1x4 Layout)');
% % xlabel('x (m)');
% % ylabel('y (m)');
% 
% % Repeat the phase screen to visualize without periodicity
% repeated_phase_screen1 = repmat(real(phz2), 4,4);
% 
% % Plot the repeated phase screen
% figure;
% imagesc(repeated_phase_screen1);
% colormap("jet");           % Apply the 'jet' colormap
% colorbar; 
% title('Repeated Randomized Phase Screen rand (1x4 Layout)');
% xlabel('x (m)');
% ylabel('y (m)');
% 
% % repeated_phase_screen3 = repmat(real(phz3), 4,4);
% % 
% % % Plot the repeated phase screen
% % figure;
% % imagesc(repeated_phase_screen3);
% % colormap("jet");           % Apply the 'jet' colormap
% % colorbar; 
% % title('Repeated Randomized Phase Screen sh (1x4 Layout)');
% % xlabel('x (m)');
% % ylabel('y (m)');
% Define different values of M (grid sizes)
M_values = [256, 512, 1024, 2048, 4096];

% Number of runs for averaging
num_runs = 5;

% Initialize arrays to store processing times for both phase screens
processing_times_phz = zeros(length(M_values), num_runs);
processing_times_phz2 = zeros(length(M_values), num_runs);

% Loop through each value of M
for i = 1:length(M_values)
    M = M_values(i);
    
    % Run multiple times for averaging PhaseScreen Traditional
    for j = 1:num_runs
        tic;
        phz = phasescreen(M, delta, L0, l0, wavelength, Cn, delta_z);
        processing_times_phz(i, j) = toc;
    end
    
    % Run multiple times for averaging PhaseScreen Randomised
    for j = 1:num_runs
        tic;
        phz2 = phasescreen2(M, delta, L0, l0, wavelength, Cn, delta_x, delta_y, delta_z);
        processing_times_phz2(i, j) = toc;
    end
end

% Calculate the average processing time for both methods
avg_processing_times_phz = mean(processing_times_phz, 2);
avg_processing_times_phz2 = mean(processing_times_phz2, 2);

% Calculate percentage difference between the average processing times
percentage_diff = 100 * (avg_processing_times_phz2 - avg_processing_times_phz) ./ avg_processing_times_phz;

% Plot the average processing times for both phase screens on the same graph
figure;
subplot(2, 1, 1); % First subplot for average processing times
plot(M_values, avg_processing_times_phz, '-o', 'DisplayName', 'PhaseScreen Traditional');
hold on;
plot(M_values, avg_processing_times_phz2, '-s', 'DisplayName', 'PhaseScreen Randomised');
xlabel('Grid Size (M)');
ylabel('Average Processing Time (seconds)');
title('Average Processing Time vs Grid Size (M)');
legend('show');
grid on;

% Plot the percentage difference in a second subplot
subplot(2, 1, 2); % Second subplot for percentage difference
plot(M_values, percentage_diff, '-d', 'DisplayName', 'Percentage Difference');
xlabel('Grid Size (M)');
ylabel('Percentage Difference (%)');
title('Percentage Difference between PhaseScreen Methods');
grid on;
