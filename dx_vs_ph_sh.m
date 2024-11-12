% Set up a phase matrix `ph`, mask `mask`, and spacing `delta`.
M = 512; % Example matrix size
delta = 1; % Example point spacing
ph1 = phasescreen(M, delta, 1000, 0.01, 5e-7,1e-7, 1);
ph2 = phasescreen2(M, delta, 1000, 0.01, 5e-7,1e-7, delta, delta, 1);
ph = phasescreen2_sh(M, delta, 1000, 0.01, 5e-7,1e-7, delta, delta, 1, 1); % Example phase matrix with random values
mask = ones(M, M); % Example mask (uniform region)

% Compute structure function `D` using the provided function
D = str_fcn2_ft(abs(ph), mask, delta);
D1 = str_fcn2_ft(abs(ph1), mask, delta);
D2 =  str_fcn2_ft(abs(ph2), mask, delta);

% Extract diagonal elements of D, D1, and D2
D_diag = diag(D);
D1_diag = diag(D1);
D2_diag = diag(D2);

% Create x-axis for plotting, limited to M/2
x_diag = (0:floor(M/2)-1) * delta; % Adjusted for indexing along the diagonal up to M/2

% Plot the diagonal elements up to M/2
figure;
hold on;
plot(x_diag, D_diag(1:floor(M/2)), 'DisplayName', 'D');
plot(x_diag, D1_diag(1:floor(M/2)), 'DisplayName', 'D1');
plot(x_diag, D2_diag(1:floor(M/2)), 'DisplayName', 'D2');
hold off;
title('Diagonal Elements of Structure Function D (up to M/2)');
xlabel('Diagonal Index (m)');
ylabel('D(x) Diagonal Values');
legend;
grid on;