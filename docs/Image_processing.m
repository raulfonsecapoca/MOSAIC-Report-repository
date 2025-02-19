%% Part Test
clear

I = imread('avenidaPaulista.jpg');
Obj = rgb2gray(I);
%imwrite(Obj, "BW.png");
% 1024 x 768 Sensor size
X = -512:1:511;
Y = -384:1:383;
[XX, YY] = meshgrid(X,Y);
dx = 12e-6;
dy = dx;
R_aperure = 4.061 * 1e-3;

% Create the PSF as a circular aperture
P = (XX * dx).^2 + (YY * dy).^2 < R_aperure^2;


% Compute the Point Spread Function (PSF) in the frequency domain
beta = 48;
chi = 0.5;
alpha = 17.3;

phase_mask_exp = beta * (XX.*dx) .* exp(chi * (XX.*dx).^2) + beta * (YY.*dy) .* exp(chi * (YY.*dy).^2);
phase_mask_cubic = alpha * (XX.*dx).^3 + alpha * (YY.*dx).^3;
lambda = 10e-6; % 8 to 12 microns
f = 14e-3;
d_i = 14e-3;


% Desired Signal-to-Noise Ratio (SNR)
SNR = 15;


% Compute the power spectral densities
S_OO = abs(fft2(Obj)).^2;  % Power spectral density of original image
S_nn = S_OO / (10^(SNR / 10));


% Get the size of the image
[rows, cols] = size(Obj);

% Generate a noise matrix with Gaussian distribution (mean 0, std 1)
noise_matrix = randn(rows, cols);

% Scale the noise to achieve the desired noise power
noise_matrix = noise_matrix .* sqrt(S_nn);

sum_fftPSF_exp = 0;
sum_fftPSF_cubic = 0;
sum_fftPSF_square_exp = 0;
sum_fftPSF_square_cubic = 0;
sum_MSE_cubic = 0;


%%%%% reading psf data%%%%%%%%%%%%%%%%%%%%%%



file_names_cubic = ["PSF_cubic_0(3)_440(9).xlsx","PSF_cubic_0(5)_438(1).xlsx", ...
    "PSF_cubic_1_436(1).xlsx", "PSF_cubic_2m.xlsx", "PSF_cubic_4m.xlsx", "PSF_cubic_6m.xlsx", ...
    "PSF_cubic_500_434(2).xlsx", "PSF_1E6.xlsx", "PSF_cubic.xlsx"];
file_names_exp = ["PSF_exp_1E6.xlsx", "PSF_exp_7.xlsx", "PSF_exp_1E6.xlsx", "PSF_exp_1E6.xlsx", ...
    "PSF_exp_1E6.xlsx", "PSF_exp_1E6.xlsx", "PSF_exp_1E6.xlsx", "PSF_exp_1E6.xlsx", "PSF_exp_1E6.xlsx"];
size_vector = [440.9, 438.1, 436.1, 435.1, 434.7, 434.5, 434.2, 434.5, 434.5];

n = length(file_names_cubic);

%real_PSF_cubic;

filename_cubic = "PSF_1E6.xlsx";
PSF_origin_cubic = readmatrix(filename_cubic);

v = linspace(-434.5, 434.5, 512) * 1e-6;
u = v;
[VV, UU] = meshgrid(v, u);

PSF_origin_cubic = interp2(VV, UU, PSF_origin_cubic, XX * dx, YY * dy, 'linear', 0);
PSF_origin_cubic = PSF_origin_cubic / sum(PSF_origin_cubic(:));
OTF_origin = fft2(PSF_origin_cubic);


for k = 1:n

    filename_cubic = file_names_cubic(k);
    real_PSF_cubic = readmatrix(filename_cubic);

    filename_exp = file_names_exp(k);
    real_PSF_exp = readmatrix(filename_exp);

    v_cubic = linspace(-size_vector(k), size_vector(k), 512) * 1e-6;
    u_cubic = v_cubic;
    [VV_cubic, UU_cubic] = meshgrid(v_cubic, u_cubic);

    v_exp = linspace(-436.9, 436.9, 512) * 1e-6;
    u_exp = v_exp;
    [VV_exp, UU_exp] = meshgrid(v_exp, u_exp);

    real_PSF_cubic = interp2(VV_cubic, UU_cubic, real_PSF_cubic, XX * dx, YY * dy, 'linear', 0);
    real_PSF_cubic = real_PSF_cubic / sum(real_PSF_cubic(:));

    real_PSF_exp = interp2(VV_exp, UU_exp, real_PSF_exp, XX * dx, YY * dy, 'linear', 0);
    real_PSF_exp = real_PSF_exp / sum(real_PSF_exp(:));

    PSF_cubic = real_PSF_cubic;
    fftPSF_cubic = fft2(PSF_cubic);
    fftPSF_cubic_centered = calibration(OTF_origin, fftPSF_cubic, XX * dx, YY * dy);

    PSF_exp = real_PSF_exp;
    fftPSF_exp = fft2(PSF_exp);
    


    
    sum_fftPSF_cubic = sum_fftPSF_cubic + conj(fftPSF_cubic_centered);
    sum_fftPSF_square_cubic = sum_fftPSF_square_cubic + abs(fftPSF_cubic_centered).^2;

    sum_fftPSF_exp = sum_fftPSF_exp + conj(fftPSF_exp);
    sum_fftPSF_square_exp = sum_fftPSF_square_exp + abs(fftPSF_exp).^2;
 
    
end




% Compute the average Wiener filter

filter_cubic = (1/n) * sum_fftPSF_cubic ./ ((1/n) * sum_fftPSF_square_cubic + S_nn./S_OO);
filter_exp = (1/n) * sum_fftPSF_exp ./ (( 1/n) * sum_fftPSF_square_exp + S_nn./S_OO);

%%









%filter_cubic = ( conj(fft2(real_PSF_cubic)) .* S_OO ) ./ ( abs(fft2(real_PSF_cubic)).^2 .* S_OO + S_nn );
%filter_cubic = fftshift(filter_cubic);

%%


%real_PSF_exp;

filename_exp = "PSF_exp_1E6.xlsx";
real_PSF_exp = readmatrix(filename_exp);

v = linspace(-436.9, 436.9, 512) * 1e-6;
u = v;
[VV, UU] = meshgrid(v, u);

real_PSF_exp = interp2(VV, UU, real_PSF_exp, XX * dx, YY * dy, 'linear', 0);
real_PSF_exp = real_PSF_exp / sum(real_PSF_exp(:));


% Add noise to the image (no convolution at this point)

noisy_image_fft_cubic = fft2(Obj).* fft2(PSF_origin_cubic) + noise_matrix;
noisy_image_cubic = ifft2(noisy_image_fft_cubic);
noisy_image_cubic = abs(fftshift(noisy_image_cubic));

noisy_image_fft_exp = fft2(Obj).* fft2(real_PSF_exp) + noise_matrix;
noisy_image_exp = ifft2(noisy_image_fft_exp);
noisy_image_exp = abs(fftshift(noisy_image_exp));




% Apply the Wiener filter by multiplying in the frequency domain
recovered_fft_cubic = noisy_image_fft_cubic .* filter_cubic;
recovered_fft_exp = noisy_image_fft_exp .* filter_exp;

% Get the recovered image by inverse FFT
recovered_Image_cubic = abs(ifft2(recovered_fft_cubic));
recovered_Image_exp = abs(ifft2(recovered_fft_exp));


% Visualization part:

% Normalize the noisy and recovered images for proper display
noisy_image_cubic = noisy_image_cubic - min(noisy_image_cubic(:));  % Shift so that the minimum is zero
noisy_image_cubic = noisy_image_cubic / max(noisy_image_cubic(:)) * 255;  % Scale to the range [0, 255]

recovered_Image_cubic = recovered_Image_cubic - min(recovered_Image_cubic(:));  % Normalize recovered image
recovered_Image_cubic = recovered_Image_cubic / max(recovered_Image_cubic(:)) * 255;  % Scale to the range [0, 255]


noisy_image_exp = noisy_image_exp - min(noisy_image_exp(:));  % Shift so that the minimum is zero
noisy_image_exp = noisy_image_exp / max(noisy_image_exp(:)) * 255;  % Scale to the range [0, 255]

recovered_Image_exp = recovered_Image_exp - min(recovered_Image_exp(:));  % Normalize recovered image
recovered_Image_exp = recovered_Image_exp / max(recovered_Image_exp(:)) * 255;  % Scale to the range [0, 255]

MSE_recovered = immse(recovered_Image_cubic, double(Obj));
MSE_blurred = immse(noisy_image_cubic, double(Obj));

% Display the original image
figure(1);
subplot(2, 3, 1);
imshow(uint8(Obj)), title('Original Image');
% Display the original image
subplot(2, 3, 4);
imshow(uint8(Obj)), title('Original Image');

% Display the noisy image
subplot(2, 3, 2);
imshow(uint8(noisy_image_cubic)), title('Blurred Noisy Image - Cubic phase mask');
% Display the noisy image
subplot(2, 3, 5);
imshow(uint8(noisy_image_exp)), title('Blurred Noisy Image - Exp phase mask');

% Display the recovered image
subplot(2, 3, 3);
imshow(uint8(recovered_Image_cubic)), title('Recovered Image - Cubic phase mask');
% Display the recovered image
subplot(2, 3, 6);
imshow(uint8(recovered_Image_exp)), title('Recovered Image - Exp phase mask');



% recenter PSF
function OTF = calibration(OTF_origin,OTF_calibration,u,v)
    p = struct();
    p.recalage = true; % activate calibration
    p.posFoc = 1; % Set the first PSF to be the standard of calibration

    FTO2D(:, :, 1, 1) = OTF_origin;
    FTO2D(:, :, 2, 1) = OTF_calibration;

    [FTO2D_rec,liste_dec] = recaleOTFs(FTO2D, p, u,v);

    OTF = FTO2D_rec( :, :,2);
end