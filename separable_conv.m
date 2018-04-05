% Script to test separable 2D filtering
% Brendan Lynskey 2018


clear all;
close all;
figure_num = 1;



scale_fact      = 2;
target_atten_dB = 70;

% Load test image (what else!?)
lenna = imread('lenna.bmp');

figure(figure_num); figure_num = figure_num + 1;
imshow(lenna)
title('Original Lenna');


% Resize using Octave function
figure(figure_num); figure_num = figure_num + 1;
imshow(imresize(lenna./max(lenna), scale_fact))
title('Lenna, scaled by Octave funciton');

% Design a filter
mag_sband = power(10, (target_atten_dB/-20));
f = [0, 0.9/scale_fact, 1.2/scale_fact, 1];
m = [1 1 mag_sband mag_sband];
filt_interp = firls(32, f, m);
filt_interp = filt_interp ./ sum(filt_interp); % Normalise imp-response

figure(figure_num); figure_num = figure_num + 1;
periodogram(filt_interp)



% Perform vertical convolution
mat_v_filt = [];
for idx_col = 1:size(lenna)(2)
    mat_v_filt = [mat_v_filt, conv(lenna(:, idx_col), filt_interp)];
end

% Perform horizontal convolution
mat_vh_filt = [];
for idx_row = 1:size(mat_v_filt)(1)
    mat_vh_filt = [mat_vh_filt; conv(mat_v_filt(idx_row, :)', filt_interp)'];
end

% Display the result, warts and all!
figure(figure_num); figure_num = figure_num + 1;
imshow(mat_vh_filt)
title('Lenna, scaled by my funciton');


 
 
 