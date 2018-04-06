% Polyphase Lanczos-windowed sinc interpolation filter generator
% Brendan Lynskey 2018
%


clear all
close all
figure_num = 1;

% Windowed sinc filter:
% Sinc good as a filter as symmetric, and linear phase preserves edges
% Octave sinc() function is the normalised sinc: sin(pi*x)/(pi*x)
%
% Lanczos window:
%    Single main sinc lobe
%    Zero @ lower & upper edges of main lobe, and beyond
%
% Make simple odd-order FIR, which captures central lobe & extra lobes on *both* sides.
% Odd order filter co-sites (some) output samples with the input samples


% Load filter spec
spec_filt

% Load Rec601 filter specs
spec_rec601


% PART 1
% Generate 1D windowed sinc filter (Lanczos window), over-sampled for upscaling
%
% 
% Find highest filter Fcutoff for desired attenuation from f_trans_end to PI
%     sinc_pband_scale: sinc argument scale-factor. 
%                       changes Fcutoff of sinc-filter (1 = all-pass)
%
% Filter should, after up-sampling, attenuate images sufficiently
% for given lanczos_order and intrp_ratio

lanczos_order = fir_ord_on2;

for sinc_pband_scale = 1.0:-0.001: 0.01

  % Generate normalised sinc filter
  %
  
  idx_sinc        = (-1*lanczos_order) : 1/intrp_ratio : lanczos_order;

  sinc_func       = sinc(idx_sinc .* sinc_pband_scale);
  wind_func       = sinc(idx_sinc ./ lanczos_order);

  % Windowing (element-wise mult) produces filter Impulse Response
  wsinc_func      = sinc_func .* wind_func;
    
  [Pxx, W]        = periodogram(wsinc_func'); % Calc PSD for analysis

  num_bins_2pi    = 2*(length(Pxx)-1);
  bin_trans_end   = round(num_bins_2pi*f_trans_end/intrp_ratio);
  min_atten       = max(Pxx(bin_trans_end:end));
  rel_atten_dB    = 10* log10(Pxx(1)) - 10*log10(min_atten);
  
  if (rel_atten_dB >= atten_trans_end_dB)
  
      % The spec is met if difference between DC attenuation, and
      % attenuation at f_trans_end/intrp_ratio is > required attenuation 
      printf('Check relative atten at new Fnyq * %f\n', (f_trans_end/intrp_ratio))

      printf('Required relative attenuation (in dB):\n')
      atten_trans_end_dB
      
      printf('Actual relative attenuation (min, in dB):\n')
      rel_atten_dB
  
      break
  end

end

 
% Plot sinc, window, windowed-sinc and its PSD
figure(figure_num); figure_num = figure_num + 1;
subplot(2,2,1);
plot(idx_sinc, sinc_func);
title('sinc');
subplot(2,2,2);
plot(idx_sinc, wind_func);
title('window');
subplot(2,2,3)
plot(idx_sinc, wsinc_func);
title('windowed sinc');
subplot(2,2,4);
periodogram(wsinc_func')



% PART 2
% Generate polyphase coefficients for Lanczos anti-imaging filter
num_taps        = (2*lanczos_order) + 1;
num_phases      = intrp_ratio;

fir_poly        = zeros(num_phases, num_taps);

% Zero-pad impulse response so can make equal length sub-filters
zpad_wsinc_func = [wsinc_func, zeros(1, num_phases-1)];

for idx_phs = 1:num_phases
    % Extract sub-filter for this phase
    fir_poly(idx_phs, :) = num_phases * downsample(zpad_wsinc_func, intrp_ratio, idx_phs-1);
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
fid = fopen('poly_lanczos_wsinc_coeffs.txt', 'w');
fprintf(fid, '%d\n', num_phases);
fprintf(fid, '%d\n', num_taps);
for idx_phase = 1:size(fir_poly)(1)
    for idx_tap = 1:size(fir_poly)(2)
        fprintf(fid, '%f\n', fir_poly(idx_phase, idx_tap));
    end
end
fclose(fid);


% PART 3: test
test_filter
