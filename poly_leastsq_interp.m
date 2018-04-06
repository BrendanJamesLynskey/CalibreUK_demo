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

fir_imp_resp    = ls_filt';
filename_coeffs = 'poly_leastsq_coeffs.txt';
gen_polyphase


% PART 3: test
test_filter
