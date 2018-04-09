% Script to test separable 2D filtering
% Brendan Lynskey 2018


%clear all;
%close all;



% If no LFP defined, design one here
if (exist("fir_imp_resp", "var") == 0)

    printf('No filter defined, so generating new LPF\n')
    figure_num = 1;

    % Load filter spec
    spec_filt

    firls_order     = 2 * fir_ord_on2 * intrp_ratio;
    mag_sband       = power(10, (target_atten_dB/-20));
    
    f               = [0, 0.9/intrp_ratio, 1.2/intrp_ratio, 1];
    m               = [1 1 mag_sband mag_sband];
    
    fir_imp_resp    = firls(firls_order, f, m);
    fir_imp_resp    = fir_imp_resp ./ sum(fir_imp_resp); % Normalise imp-resp
bjl
    figure(figure_num); figure_num = figure_num + 1;
    periodogram(fir_imp_resp)

end

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
    mat_v_filt = [mat_v_filt, conv(image_double_usrc(:, idx_col), fir_imp_resp)];
end

% Perform horizontal convolution
mat_vh_filt = [];
for idx_row = 1:size(mat_v_filt)(1)
    mat_vh_filt = [mat_vh_filt; conv(mat_v_filt(idx_row, :)', fir_imp_resp)'];
end

% Convert the output to uint8
mat_vh_filt_uint8 = cast(mat_vh_filt .* power(2, pxl_depth), 'uint8');

% Display the result, warts and all!
figure(figure_num); figure_num = figure_num + 1;
imshow(mat_vh_filt_uint8)
title(sprintf('%s%s', 'Image, scaled by ', filter_name));


 