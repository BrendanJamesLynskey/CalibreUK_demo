% Figure sequencing script
% Brendan Lynskey 2018




% Following code writes current figure as PDF
% Octave bug means have to do this on latest figure
op_dir        = 'figs_poly_lanczos_wsinc_interp_4phs_7taps/'

%Use these!
%intrp_ratio
%fir_ord_on2

print_fig_num = 1
op_fname      = sprintf('%sfigure_%d.pdf', op_dir, figure);

if ((figure_num - 1) == print_fig_num)
    print(op_fname);
end

bjl


% Increment the figure number and create a new figure

figure(figure_num);
figure_num = figure_num + 1;