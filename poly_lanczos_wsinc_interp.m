% Polyphase Lanczos-windowed sinc interpolation filter generator
% Brendan Lynskey 2018
%

% Design 1D interpolation filter for SD video (Rec 601)
%
% Provides a few options fo rspecification of desired filter charachteristics
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
%
% TO DO: use firls() to find hand-optimised filter, and compare performance

clear all
close all
figure_num = 1;

% Specify filter
%
%    lanczos_order:    number of lobes each side of the main-lobe.
%                      higher order yields higher performance, at greater cost
%    intrp_ratio:      desired interpolation (scaling) factor
%                      equals number of sub-filter phases
%    signed_coeff_wid: bit-width of FIR coefficients (2's complement)
%                      must represent +/-1, so 2 MS bits to left of binary point
%
lanczos_order      = 3;
intrp_ratio        = 4;
%signed_coeff_wid   = 9; % Altera multipliers used efficiently with 9b IPs
signed_coeff_wid   = 18; % Altera multipliers used fully with 18b IPs




% PART 1
% Design a filter which, after up-sampling, attenuates images sufficiently
% for given lanczos_order and intrp_ratio
%
% Specify Rec601 presampling filter mask
f_samp_orig        = 27000000;
f_presamp_40dB     =  8000000;
% Anti-imaging filter transition band specified @ lower edge of 1st image 
f_trans_beg        = (f_samp_orig * 0.25) / f_samp_orig;
f_trans_end        = (f_samp_orig - f_presamp_40dB) / f_samp_orig;

%
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



% PART 2
% Generate 1D windowed sinc filter (Lanczos window), over-sampled for upscaling
% 
% Find highest filter Fcutoff for desired attenuation from f_trans_end to PI
%     sinc_pband_scale: sinc argument scale-factor. 
%                       changes Fcutoff of sinc-filter (1 = all-pass)
%
%
for sinc_pband_scale = 1.0:-0.001: 0.01

  %printf("Check relative atten at sinc_pband_scale = %f\n", sinc_pband_scale)

  % Generate normalised sinc filter. Sinc good as filter must be symmetric,
  % as linear phase preserves edges in images
  % Then generate Lanczos window:
  %    Single main sinc lobe, scaled_taps
  %    Zero at lower and upper extremes
  idx_sinc        = (-1*lanczos_order) : 1/intrp_ratio : lanczos_order;

  sinc_func       = sinc(idx_sinc .* sinc_pband_scale);
  wind_func       = sinc(idx_sinc ./ lanczos_order);

  % Windowing i.e. element-wise mult produces prototype filter
  wsinc_func      = sinc_func .* wind_func;
    
  [Pxx, W]        = periodogram(wsinc_func'); % Calc PSD for analysis

  num_bins_2pi    = 2*(length(Pxx)-1);
  bin_trans_end   = round(num_bins_2pi*f_trans_end/intrp_ratio);
  min_atten       = max(Pxx(bin_trans_end:end));
  rel_atten_dB    = 10* log10(Pxx(1)) - 10*log10(min_atten);
  
  if (rel_atten_dB >= atten_trans_end_dB)
  
      % The spec is met if difference between DC attenuation, and
      % attenuation at f_trans_end/intrp_ratio is > required attenuation 
      printf("Check relative atten at new Fnyq * %f\n", (f_trans_end/intrp_ratio))

      printf("Required relative attenuation (in dB):\n")
      atten_trans_end_dB
      
      printf("Actual relative attenuation (min, in dB):\n")
      rel_atten_dB
  
      break
  end

end

 
% Plot sinc, window, windowed-sinc and its PSD
figure(figure_num); figure_num = figure_num + 1;
subplot(2,2,1);
plot(idx_sinc, sinc_func);
title("sinc");
subplot(2,2,2);
plot(idx_sinc, wind_func);
title("window");
subplot(2,2,3)
plot(idx_sinc, wsinc_func);
title("windowed sinc");
subplot(2,2,4);
periodogram(wsinc_func')



% PART 3
% Generate polyphase coefficients for Lanczos anti-imaging filter
% Make simple odd-order FIR, which captures central lobe,
% and the extra lobes on *both* sides.
% Odd order filter co-sites output with the input samples
num_taps        = (2*lanczos_order) + 1;
num_phases      = intrp_ratio;

fir_poly        = zeros(num_phases, num_taps);

% Generate impulse response for each of the polyphase filters
% NB: Octave sinc() function is the normalised sinc: sin(pi*x)/(pi*x)
% Integer arguments for sinc return 0 except at x=0, as each lobe occupies an
% interval of 1 (assuming passband is 100% of Nyquist band)

% Could use function downsample() here

for idx_phs = num_phases:-1:1
    for idx_tap = 1:num_taps
        % Choose x co-ordinate for sub-filter in question
        % Step over x-axis with increments of 1 (single lobe)
        % Each subsequent phase starts a fraction of a lobe to the left,
        % that's because in convoution, impulse-response is LR flipped
        phs_offset                 = (idx_phs - num_phases - 1) / num_phases;
        x_coord                    = (-1 * lanczos_order) + (idx_tap-1) + phs_offset;
        sinc_val                   = sinc(x_coord .* sinc_pband_scale);
        % Lanczos window uses a dilated central lobe to weight the sinc function,
        % weights decreasing as move away from x=0
        wind_val                   = sinc(x_coord ./ lanczos_order);
        fir_poly(idx_phs, idx_tap) = sinc_val * wind_val;
    end
    % Normalise mag-response
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






% PART 4
% Demonstrate scaling with square-wave input

ip_len = 32;
ip_sig = [ones(1,ip_len/4), zeros(1,ip_len/4), ones(1,ip_len/4), zeros(1,ip_len/4)];


% Calculate output of each phase
% Plot input signal, and results of simple filtering (sanity check)

figure(figure_num); figure_num = figure_num + 1;

fig_xaxis        = [0:(ip_len-1)];

subplot(((num_phases/2)+1),2,1);
plot(fig_xaxis, ip_sig);
axis([0, (ip_len-1), -0.5, 1.5]);

phase_op_len     = ip_len + (num_taps-1);
op_phase         = zeros(num_phases, phase_op_len);

for phase = 1:num_phases
    op_phase(phase, :) = conv(ip_sig, fir_poly(phase, :));
    subplot(((num_phases/2)+1),2,phase+2);
    plot(fig_xaxis, op_phase(phase, 1:ip_len));
    axis([0, (ip_len-1), -0.5, 1.5]);
end


% Create upscaled output, by commutating between phases
op_upsc          = zeros(1, num_phases * ip_len);

idx_phase = 1;
idx_samp  = 1;
for idx_upsc = 1: length(op_upsc)
    op_upsc(idx_upsc) = op_phase(idx_phase, idx_samp);
    if (idx_phase == num_phases)
        idx_phase         = 1;
        idx_samp          = idx_samp + 1;
    else
        idx_phase         = idx_phase + 1;
    end
end

% Plot the input and output signals. NB group-delay and x-axis scale
figure(figure_num); figure_num = figure_num + 1;

subplot(2,1,1);
plt_xaxis        = [0:(ip_len-1)];
plot(plt_xaxis, ip_sig);
axis([0, (ip_len-1), -0.5, 1.5]);

subplot(2,1,2);
plt_xaxis        = [0:(length(op_upsc)-1)];
plot(plt_xaxis, op_upsc);
