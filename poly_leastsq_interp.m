% Least squares approximaiton interpolation filter generator
% Brendan Lynskey 2018
%


clear all
close all
figure_num = 1;

% Least-squares filter:
%  This is just an FIR filter designed via the Octave firls() function
%  From the help:
%     FIR filter design using least squares method.  Returns a length N+1
%     linear phase filter such that the integral of the weighted mean
%     squared error in the specified bands is minimized.
% This is good as linear phase preserves edges


% Load filter spec
spec_filt

% Load Rec601 filter specs
spec_rec601


% PART 1
% Generate 1D least-squares approx filter, over-sampled for upscaling


BW_old_Fnyq     = f_presamp_40dB/(f_samp_orig/2);
BW_new_Fnyq     = BW_old_Fnyq/intrp_ratio;


% Split Nyquist band into bands, forming freq and mag vectors
f_vect    = [];
m_vect    = [];
mag_sband = power(10, (atten_trans_end_dB/-20));

for idx_phs = 1 : (intrp_ratio/2+1)
    mid_f   = (idx_phs-1) * 2/intrp_ratio;
    lower_f = mid_f - BW_new_Fnyq;
    upper_f = mid_f + BW_new_Fnyq;
    f_vect = [f_vect, max(0, lower_f), min(1, upper_f)];
    if (idx_phs == 1)
        m_vect = [m_vect, 1, 1];
    else
        m_vect = [m_vect, mag_sband, mag_sband];
    end
end

% TODO: could play with relative-weighting for each band
%w_vect=[...];


ls_filt = firls(intrp_ratio*fir_ord_on2*2,f_vect,m_vect);

% Use following command to ensure that impulse response is symmetric
%ls_filt = (ls_filt+fliplr(ls_filt))/2;
 
% Plot filter and its PSD
figure(figure_num); figure_num = figure_num + 1;
subplot(2,1,1);
stem(ls_filt);
title('Least squares approx filter');
subplot(2,1,2);
periodogram(ls_filt')




% PART 2
% Generate polyphase coefficients for LS approx anti-imaging filter
num_taps        = (fir_ord_on2*2) + 1;
num_phases      = intrp_ratio;

fir_poly        = zeros(num_phases, num_taps);

% Zero-pad impulse response so can make equal length sub-filters
zpad_ls_filt    = [ls_filt', zeros(1, num_phases-1)];

for idx_phs = 1:num_phases
    sub_filter              = downsample(zpad_ls_filt, num_phases, idx_phs-1)';
    fir_poly(idx_phs, :)    = num_phases .* sub_filter;
    % Normalise imp-response
    sum_mag                 = sum(fir_poly(idx_phs, :));
    fir_poly(idx_phs, :)   /= sum_mag;
    % Round coefficients to nearest fixed-point value
    % For range of +/-1, need 2 integer bits
    scale_fact              = 2 ^ (signed_coeff_wid - 2);
    scaled_taps             = fir_poly(idx_phs, :) .* scale_fact;
    fir_poly(idx_phs, :)    = round(scaled_taps)   ./ (scale_fact);
end


% Plot PSD of each phase. Check that all have low-pass spectrum!
figure(figure_num); figure_num = figure_num + 1;

subplot_num = 1;
for phase = 1 : num_phases
    subplot(num_phases, 2, subplot_num);
    subplot_num = subplot_num + 1;
    plot(1:num_taps, fir_poly(phase, :));

    subplot(num_phases, 2, subplot_num);
    subplot_num = subplot_num + 1;
    periodogram(fir_poly(phase, :)');
end

% Save coefficients in file
fid = fopen('poly_leastsq_coeffs.txt', 'w');
fprintf(fid, '%d\n', num_phases);
fprintf(fid, '%d\n', num_taps);
for idx_phs = 1:size(fir_poly)(1)
    for idx_tap = 1:size(fir_poly)(2)
        fprintf(fid, '%f\n', fir_poly(idx_phs, idx_tap));
    end
end
fclose(fid);




% PART 4: test
test_filter
