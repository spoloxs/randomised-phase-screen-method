% Parameters
delta = 1;   % Grid spacing [m]
M = 512;     % Grid size (constant)
l0 = 1;      % Inner scale [m]
delta_x = delta;
delta_y = delta;
delta_z = delta;
wavelength = 5e-7;
Cn = 1e-7;
L0_values = [M*delta_x, 100*M*delta_x, 10000*M*delta_x, 100000*M*delta_x];  % Different outer scale values [m]

% Define parameters for phase screen generation
num_trials = 500;  % Reduced trial count for computational efficiency
central_phase_samples = zeros(num_trials, length(L0_values));
relative_phase_samples = zeros(num_trials, length(L0_values));
log_magnitude_phase_samples = zeros(num_trials, length(L0_values));

% Loop over different values of L0 and calculate phase samples
for L0_idx = 1:length(L0_values)
    L0 = L0_values(L0_idx);  % Update L0 for each iteration
    
    for trial = 1:num_trials
        % Step 1: Generate random phase screen
        phase_screen = phasescreen2(M, delta, L0, l0, wavelength, Cn, delta_x, delta_y, delta_z);
        
        % Step 2: Log key parameters
        center_phase = phase_screen(M/2, M/2);
        relative_phase = phase_screen(M/2, M/2) - mean(phase_screen(:));
        log_magnitude_phase = log(abs(center_phase));
        
        % Store results for each trial
        central_phase_samples(trial, L0_idx) = center_phase;
        relative_phase_samples(trial, L0_idx) = relative_phase;
        log_magnitude_phase_samples(trial, L0_idx) = log_magnitude_phase;
    end
end

% Plotting all histograms in subplots
figure;
set(gcf, 'Position', [100, 100, 1400, 800])
for L0_idx = 1:length(L0_values)
    L0 = L0_values(L0_idx);
    
    % Subplot for absolute phase
    subplot(3, length(L0_values), L0_idx);
    histogram(real(central_phase_samples(:, L0_idx)), 'Normalization', 'pdf');
    title(['Abs Phase L0 = ', num2str(L0)]);
    xlabel('Phase');
    ylabel('PDF');
    
    % Subplot for relative phase
    subplot(3, length(L0_values), length(L0_values) + L0_idx);
    histogram(real(relative_phase_samples(:, L0_idx)), 'Normalization', 'pdf');
    title(['Rel Phase L0 = ', num2str(L0)]);
    xlabel('Phase');
    ylabel('PDF');
    
    % Subplot for log magnitude of absolute phase
    subplot(3, length(L0_values), 2*length(L0_values) + L0_idx);
    histogram(real(log_magnitude_phase_samples(:, L0_idx)), 'Normalization', 'pdf');
    title(['Log Abs Phase L0 = ', num2str(L0)]);
    xlabel('Log Phase');
    ylabel('PDF');
end

% Adjust layout for better visibility
sgtitle('Phase Distribution for Different L0 Values');

