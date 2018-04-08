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

filter_name = 'lanczos-winsinc-filter';

% Load filter spec
spec_filt

% Load Rec601 signal spec
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
  fir_imp_resp      = sinc_func .* wind_func;
    
  [Pxx, W]        = periodogram(fir_imp_resp'); % Calc PSD for analysis

  num_bins_2pi    = 2*(length(Pxx)-1);
  bin_trans_end   = round(num_bins_2pi*f_trans_end/intrp_ratio);
  
  % Currenly check all bands above end of transition band.
  % TODO: could only check bands which can have significant image power.
  %       Probably wouldn't help as stopband lobes are wide,
  %       and they become smaller in magnitude at higher frequencies
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
      
      % Normalise imp-response magnitude
      fir_imp_resp /= sum(fir_imp_resp);
  
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
plot(idx_sinc, fir_imp_resp);
title('windowed sinc');
subplot(2,2,4);
periodogram(fir_imp_resp')



% PART 2
% Generate polyphase coefficients for Lanczos anti-imaging filter
gen_polyphase

% PART 3: test
test_1D
test_2D
