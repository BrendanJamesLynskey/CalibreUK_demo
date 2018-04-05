% Least squares approximaiton interpolation filter generator
% Brendan Lynskey 2018
%

% Design 1D interpolation filter for SD video (Rec 601)
%
% Provides a few options for specification of desired filter charachteristics
%
%
% Real video signals captured using pre-sampling anti-alias filter
% BT.601 recommended pre-sampling filter is ~halfband, as follows:
%    12dB atten @6.75MHz (@0.25 x Fs)
%    40dB atten @8.00MHz (@0.30 x Fs)
%
% Will assume all significant signal-energy contained within this mask
% Aim to make SNR no worse than this
%
% Dynamic range of signals, assuming that use full-scale:
%   8b samples: ~48dB
%  10b samples: ~60dB
%
%
% NB: not specifying passband ripple
%

clear all
close all
figure_num = 1;

% Specify filter
%
%    order_on_2:       number of taps each side of central tap
%                      higher order yields higher performance, at greater cost
%    intrp_ratio:      desired interpolation (scaling) factor
%                      equals number of sub-filter phases
%    signed_coeff_wid: bit-width of FIR coefficients (2's complement)
%                      must represent +/-1, so 2 MS bits to left of binary point
%
order_on_2         = 3;
intrp_ratio        = 4;
%signed_coeff_wid   = 9; % Altera multipliers used efficiently with 9b IPs
signed_coeff_wid   = 18; % Altera multipliers used fully with 18b IPs

% Load Rec601 filter specs
spec_rec601



% PART 1
% Generate 1D least-squares approx filter, over-sampled for upscaling
num_taps        = (2*order_on_2) + 1;
num_phases      = intrp_ratio;

BW_old_Fnyq     = f_presamp_40dB/(f_samp_orig/2);
BW_new_Fnyq     = BW_old_Fnyq/num_phases;


% Split Nyquist band into bands, forming freq and mag vectors
f_vect = [];
m_vect = [];

mag_sband = power(10, (atten_trans_end_dB/-20));

for idx_phs = 1 : (num_phases/2+1)
    mid_f   = (idx_phs-1) * 2/num_phases;
    lower_f = mid_f - (BW_new_Fnyq);
    upper_f = mid_f + (BW_new_Fnyq);   
    f_vect = [f_vect, max(0, lower_f), min(1, upper_f)];
    if (idx_phs == 1)
        m_vect = [m_vect, 1, 1];
    else
        m_vect = [m_vect, mag_sband, mag_sband];
    end
end

% TODO: could play with relative-weighting for each band
%w_vect=[...];


ls_filt = firls(num_phases*num_taps,f_vect,m_vect);

% Use this command to ensure that impulse response is symmetric
%ls_filt = (ls_filt+fliplr(ls_filt))/2;
 
% Plot filter and its PSD
figure(figure_num); figure_num = figure_num + 1;
subplot(2,1,1);
stem(ls_filt);
title("Least squares approx filter");
subplot(2,1,2);
periodogram(ls_filt')




% PART 2
% Generate polyphase coefficients for LS approx anti-imaging filter

% Zero-pad taps to power-of-two, including central impulse
num_taps        = num_taps + 1;

fir_poly        = zeros(num_phases, num_taps);

zeropad_pow2    = ceil(log2(length(ls_filt)));
ls_filt_zpadlen = pow2(zeropad_pow2) - length(ls_filt);
ls_filt_zpad    = [ls_filt', zeros(1, ls_filt_zpadlen)];

for idx_phs = 1:num_phases
    fir_poly(idx_phs, :) = downsample(ls_filt_zpad, num_phases, idx_phs-1)';    
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
fprintf(fid, "%d\n", num_phases);
fprintf(fid, "%d\n", num_taps);
for idx_phs = 1:size(fir_poly)(1)
    for idx_tap = 1:size(fir_poly)(2)
        fprintf(fid, "%f\n", fir_poly(idx_phs, idx_tap));
    end
end
fclose(fid);




% PART 4: test
test_filter
