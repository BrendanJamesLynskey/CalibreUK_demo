% Calibre UK video scaling demo
% Brendan Lynskey 2018

setup_octave

clear all
close all

en_write_plot_pdf          = 0;

gen_lanczos_winsinc_filter = 0;
gen_least_sq_approx_filter = 1;


% Generate Lanczos windowed-sinc filter
%
if (gen_lanczos_winsinc_filter == 1)
    % Setup
    keep gen_lanczos_winsinc_filter gen_least_sq_approx_filter en_write_plot_pdf;
    close all
    figure_num = 1;
    % Load filter spec
    spec_filt
    % Load Rec601 signal spec
    spec_rec601
    % Generate the filter, storing plots as PDF
    filter_name = 'lanczos-winsinc-filter';
    poly_lanczos_wsinc_interp
    % Generate polyphase coefficients for LS approx anti-imaging filter
    gen_polyphase
    % Test
    test_1D
    test_2D
end





% Generate least-squares approx filter
%
if (gen_least_sq_approx_filter == 1)
    % Setup
    keep gen_lanczos_winsinc_filter gen_least_sq_approx_filter en_write_plot_pdf;
    close all
    figure_num = 1;
    % Load filter spec
    spec_filt
    % Load Rec601 signal spec
    spec_rec601
    % Generate the filter, storing plots as PDF
    filter_name = 'leastsq-approx-filter';
    poly_leastsq_interp
    % Generate polyphase coefficients for LS approx anti-imaging filter
    gen_polyphase
    % Test
    test_1D
    test_2D
end

