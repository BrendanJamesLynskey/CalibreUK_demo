% Script to specify filter
% Brendan Lynskey 2018


% Params:
%
%    fir_ord_on2:      half num taps in symm-FIR (excluding central impulse)
%                      higher order yields higher performance, at greater cost
%    intrp_ratio:      desired interpolation (scaling) factor
%                      equals number of sub-filter phases
%    signed_coeff_wid: bit-width of FIR coefficients (2's complement)
%                      must represent +/-1, so 2 MS bits to left of binary point
%

%fir_ord_on2         = 2;   % Lanczos2 fails to meet attenuation spec.
                            % firls() produces bad filter
fir_ord_on2         = 3;    % Both filters good at this spec

%intrp_ratio         = 2;
intrp_ratio         = 4;

%signed_coeff_wid    = 9; % Altera multipliers used efficiently with 9b IPs
signed_coeff_wid    = 18; % Altera multipliers used fully with 18b IPs

% Two options wrt location of first OP-sample:
%   1) co-sited with central central IP-sample
%   2) sited 1/2 an OP-sample before  central central IP-sample
cosite_op1_at_ip    = 1;