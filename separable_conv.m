% Script to test separable 2D filtering
% Brendan Lynskey 2018


clear all;
close all;
figure_num = 1;


intrp_ratio     = 4;
target_atten_dB = 70;

% Design an LPF
firls_order     = 28;
mag_sband       = power(10, (target_atten_dB/-20));
f               = [0, 0.9/intrp_ratio, 1.2/intrp_ratio, 1];
m               = [1 1 mag_sband mag_sband];
filt_interp     = firls(firls_order, f, m);
filt_interp     = filt_interp ./ sum(filt_interp); % Normalise imp-response

figure(figure_num); figure_num = figure_num + 1;
periodogram(filt_interp)


% Load 8b greyscale test-image
pxl_depth       = 8;
image           = imread('lenna_256x256.bmp');
figure(figure_num); figure_num = figure_num + 1;
imshow(image)
title('Original image');


% Resize using Octave function
figure(figure_num); figure_num = figure_num + 1;
imshow(imresize(image, intrp_ratio))
title('Image, scaled by Octave function');

% Convert the output to doubles
image_double    = cast(image, 'double') ./ power(2, pxl_depth);

% Upsample the image
image_double_usc  = upsample(image_double,      intrp_ratio);  % Column upsample
image_double_usc  = image_double_usc  .* intrp_ratio;          % Mag scale
image_double_usrc = upsample(image_double_usc', intrp_ratio)'; % Row upsample
image_double_usrc = image_double_usrc .* intrp_ratio;          % Mag scale

% Perform vertical convolution
mat_v_filt = [];
for idx_col = 1:size(image_double_usrc)(2)
    mat_v_filt = [mat_v_filt, conv(image_double_usrc(:, idx_col), filt_interp)];
end

% Perform horizontal convolution
mat_vh_filt = [];
for idx_row = 1:size(mat_v_filt)(1)
    mat_vh_filt = [mat_vh_filt; conv(mat_v_filt(idx_row, :)', filt_interp)'];
end

% Convert the output to uint8
mat_vh_filt_uint8 = cast(mat_vh_filt .* power(2, pxl_depth), 'uint8');

% Display the result, warts and all!
figure(figure_num); figure_num = figure_num + 1;
imshow(mat_vh_filt_uint8)
title('Image, scaled by my function');


 
 
 