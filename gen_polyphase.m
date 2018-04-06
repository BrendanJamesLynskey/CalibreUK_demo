% Script to decompose an FIR impulse response into a polyphase filter
% Brendan Lynskey 2018

num_taps        = (2*fir_ord_on2) + 1;
num_phases      = intrp_ratio;

fir_poly        = zeros(num_phases, num_taps);

% Zero-pad impulse response so can make equal length sub-filters
zpad_wsinc_func = [fir_imp_resp, zeros(1, num_phases-1)];

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
fid = fopen(filename_coeffs, 'w');
fprintf(fid, '%d\n', num_phases);
fprintf(fid, '%d\n', num_taps);
for idx_phase = 1:size(fir_poly)(1)
    for idx_tap = 1:size(fir_poly)(2)
        fprintf(fid, '%f\n', fir_poly(idx_phase, idx_tap));
    end
end
fclose(fid);

