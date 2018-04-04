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

% Specify required attenuation. Build up from abse-spec
% Spec image power @end of transition band (worst-case folding here)
%
% Option 1) Spec attenuation guaranteed by Rec601 presampling filter (40dB)
atten_trans_end_dB = 40;
% Option 2) Spec SNR implied by 8b dynamic range (~48dB)
%atten_trans_end_dB = 6.02 * 8;
% Option 3) Spec SNR implied by 10b dynamic range (~60dB)
%atten_trans_end_dB = 6.02 * 10;
%
% Add extra 3dB to make folded-down power insignificant
atten_trans_end_dB = atten_trans_end_dB + 3;
% Add factor for number of images which will fold-down                          - CHECK!
atten_trans_end_dB = atten_trans_end_dB + log10(intrp_ratio-1);



% Load Rec601 specs
spec_rec601



% PART 1
% Generate 1D least-squares approx filter, over-sampled for upscaling

num_taps        = (2*order_on_2) + 1;
num_phases      = intrp_ratio;

BW_old_Fnyq     = f_presamp_40dB/(f_samp_orig/2)
BW_new_Fnyq     = BW_old_Fnyq/num_phases


% Split Nyquist band into bands, forming freq and mag vectors
f_vect = [];
m_vect = [];
for idx_phase = 1 : (num_phases+1)
    mid_f   = (idx_phase-1) * 1/num_phases
    lower_f = mid_f - (BW_new_Fnyq/2);
    upper_f = mid_f + (BW_new_Fnyq/2);    
    f_vect = [f_vect, max(0, lower_f), min(1, upper_f)];
    if (idx_phase == 1)
        m_vect = [m_vect, 1, 1];
    else
        m_vect = [m_vect, 0, 0];
    end
end


% Relative weighting for each band
%w=[1 4000000 2000000 200000 200000 100000 100000 100000 ];

ls_filt = firls(num_phases*num_taps,f_vect,m_vect);
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

fir_poly        = zeros(num_phases, num_taps);

size(ls_filt)
size(fir_poly(1, :))
myVar = downsample(ls_filt, num_phases, 0)(1:num_taps)'
size(myVar)


for idx_phase = 1:num_phases
    phase_imp_resp         = downsample(ls_filt, num_phases, idx_phase-1);
    fir_poly(idx_phase, :) = phase_imp_resp(1:num_taps)';
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
for idx_phase = 1:size(fir_poly)(1)
    for idx_tap = 1:size(fir_poly)(2)
        fprintf(fid, "%f\n", fir_poly(idx_phase, idx_tap));
    end
end
fclose(fid);




% PART 4: test
test_filter
