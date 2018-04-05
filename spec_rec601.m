% Specify Rec601 presampling filter
% Brendan Lynskey 2018

% Specify Rec601 presampling filter mask
f_samp_orig        = 27000000;
f_presamp_40dB     =  8000000;
% Anti-imaging filter transition band specified @ lower edge of 1st image 
f_trans_beg        = (f_samp_orig * 0.25) / f_samp_orig;
f_trans_end        = (f_samp_orig - f_presamp_40dB) / f_samp_orig;



% Specify required attenuation. Build up from base-spec
% Spec image power @end of transition band (worst-case folding here)
%
% Option 1) Spec attenuation guaranteed by Rec601 presampling filter (40dB)
atten_trans_end_dB = 40;

% Option 2) Spec SNR implied by 8b dynamic range (~48dB)
%atten_trans_end_dB = 6.02 * 8;

% Option 3) Spec SNR implied by 10b dynamic range (~60dB)
%atten_trans_end_dB = 6.02 * 10;





% Add extra 3dB to make folded-down power insignificant
atten_trans_end_dB = atten_trans_end_dB + 3;
% Add factor for number of images which will fold-down
atten_trans_end_dB = atten_trans_end_dB + log10(intrp_ratio);