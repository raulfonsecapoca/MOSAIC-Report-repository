clear all; % Clear all variables from memory
close all; % Close all figure windows
clc; % Clear command window

I = imread('avenidaPaulista.jpg'); % Read input image
Obj = im2gray(I); % Convert image to grayscale

% Take meter as the unit
% Compute the Point Spread Function (PSF) in the frequency domain
dx = 1e-4; % Spatial resolution
X = -512:1:511;
Y = -384:1:383;
[XX, YY] = meshgrid(X,Y);
u = XX./512;
v = YY./384;
P = u.^2 + v.^2 < 1; % Phase mask (pupil function), u and v are normalized pupil coordinates

lambda = 10 * 1e-6; % Wavelength in meters
R_aperure = 5 * 1e-3; % Aperture radius in meters
f = 14 * 1e-3; % Focal length in meters
d_start = 0.5; % Initial distance in meters
d_o = 1; % Object distance in meters
d_i = 14 * 1e-3; % Image distance in meters
SNR = 10; % Signal-to-noise ratio

S_OO = abs(fft2(Obj)).^2; % Compute the power spectrum of the object
S_nn = S_OO / (10^(SNR / 10)); % Compute noise power spectrum

load('zernike_data.mat'); % Load Zernike coefficients
phase_zernike_data = struct([]);

for idx = 1:length(zernike_data)
    c_nm = zernike_data(idx).c_nm; % Zernike coefficients
    zernike_orders = zernike_data(idx).zernike_orders; % Zernike orders
    distance = zernike_data(idx).position; % Position (distance)
    
    phase_zernike = zeros(size(u)); % Initialize phase Zernike component
    
    for i = 1:length(c_nm)
        n = zernike_orders(i, 1);
        m = zernike_orders(i, 2);
        phase_zernike = phase_zernike + c_nm(i) * zernike_polynomial(n, m, u, v);
    end
    
    phase_zernike = 2 * pi * phase_zernike / (lambda * 10^6); % Normalize phase component
    phase_zernike_data(idx).distance = distance; % Store distance
    phase_zernike_data(idx).phase_zernike = phase_zernike; % Store phase Zernike result
end

phase_zernike = phase_zernike_data(1).phase_zernike; % Select phase Zernike at distance = 1.5m

beta = 54.1345;
chi = 0.4997;

k = 1.5;
psi = pi * R_aperure^2 / lambda * (1/f - 1/k - 1/d_i); % Compute phase shift
phase_mask = pi* (beta .* (u) .* exp(chi * (u).^2) + beta .* (v) .* exp(chi * (v).^2)); % Compute phase mask
h_psi = P .* exp(1j * (- phase_zernike + phase_mask + psi * ((u).^2 + (v).^2))); % Apply phase corrections

% Display phase mask
figure;
imagesc(mod(phase_mask,2*pi));
colormap("gray");

% 3D plot of phase mask
figure;
surf(u,v,phase_mask)
shading interp;
colormap("gray")
colorbar;
xlabel('u');
ylabel('v');
title('phase_mask');

% 3D plot of h_psi phase component
figure;
surf(u,v,angle(h_psi));
shading interp;
colormap("gray")
colorbar;
xlabel('u');
ylabel('v');
title('h_psi');

% Display phase angle of h_psi
figure;
imagesc(angle(h_psi));

% Compute and normalize the PSF
PSF = abs(fft2(h_psi)).^2;
PSF = PSF / sum(PSF(:));

% Compute the Modulation Transfer Function (MTF)
OTF = fft2(PSF);
MTF = abs(OTF);
MTF = MTF / max(MTF(:));
MTF = fftshift(MTF);

% Display MTF
figure;
imagesc(MTF);
colorbar;
title('MTF of PSF');

% Extract and plot MTF profile along the center line
[nx, ~] = size(MTF);
mid_x = floor(nx / 2);
MTF_profile = MTF(mid_x, :);

figure;
plot(MTF_profile, 'b-', 'LineWidth', 1.5); % Plot MTF
legend('MTF', 'MTF\_1');
title('Comparison of MTF Profiles (Center Horizontal Line)');
xlabel('Frequency Index');
ylabel('MTF Amplitude');
grid on;

% Compute the Fourier transform of the PSF with object size adjustment
PSF = fft2(h_psi, size(Obj, 1), size(Obj, 2));
PSF_fft = fftshift(PSF);

S_OO = abs(fft2(Obj)).^2;
S_nn = S_OO / (10^(SNR / 10));
OTF = fftshift(fft2(PSF));
MTF = abs(OTF);
MTF = fftshift(MTF / max(MTF(:)));

% Function to compute Zernike polynomials
function Z = zernike_polynomial(n, m, u, v)
    [theta, rho] = cart2pol(u, v); % Convert to polar coordinates
    rho(rho > 1) = 0; % Limit to unit disk

    R = 0;
    for k = 0:floor((n - abs(m)) / 2)
        R = R + ((-1)^k * factorial(n - k)) / ...
                (factorial(k) * factorial((n + abs(m))/2 - k) * factorial((n - abs(m))/2 - k)) * rho.^(n - 2*k);
    end

    if m >= 0
        Z = R .* cos(m * theta); % Compute Zernike function with cosine term
    else
        Z = R .* sin(abs(m) * theta); % Compute Zernike function with sine term
    end
end