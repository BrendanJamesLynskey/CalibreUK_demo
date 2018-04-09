% Interpolation filter test script
% Brendan Lynskey 2018



% Demonstrate 1D scaling with square-wave input

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

