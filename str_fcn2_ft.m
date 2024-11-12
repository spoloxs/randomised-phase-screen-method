function D = str_fcn2_ft(ph, mask, delta)
% STR_FCN2_FT computes the structure function using Fourier transforms.
% Inputs:
%   ph    - The phase matrix (2D array).
%   mask  - The mask matrix (2D array) to isolate specific regions.
%   delta - The spacing between points in the phase matrix.
% Output:
%   D     - The computed structure function.

% Get the size of the input matrix
N = size(ph, 1);

% Apply the mask to the phase matrix
ph = ph .* mask;

% Compute Fourier transforms of the phase and related terms
P = fft2(ph) * delta^2;          % Fourier transform of the phase matrix, scaled by delta^2
S = fft2(ph.^2) * delta^2;       % Fourier transform of the squared phase matrix, scaled by delta^2
W = fft2(mask) * delta^2;        % Fourier transform of the mask, scaled by delta^2

% Frequency resolution
delta_f = 1 / (N * delta);

% Compute the inverse Fourier transform of the squared mask transform
w2 = ifft2(W .* conj(W)) * (N * delta_f)^2; % w2 normalization using the inverse FT, scaled by (N * delta_f)^2

% Compute the structure function using inverse Fourier transforms
D = 2 * ifft2(real(S .* conj(W)) - abs(P).^2) * (N * delta_f)^2 ./ w2 .* mask; % Scaled by (N * delta_f)^2
end