% Specify Rec601 presampling filter
% Brendan Lynskey 2018

% Specify Rec601 presampling filter mask
f_samp_orig        = 27000000;
f_presamp_40dB     =  8000000;
% Anti-imaging filter transition band specified @ lower edge of 1st image 
f_trans_beg        = (f_samp_orig * 0.25) / f_samp_orig;
f_trans_end        = (f_samp_orig - f_presamp_40dB) / f_samp_orig;
