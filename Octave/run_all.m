% Calibre UK video scaling demo
pkg load signal
pkg load image

clear all
close all

enable_gen_lanczos_winsinc_filter = 0;
enable_gen_least_sq_approx_filter = 1;


% Generate Lanczos windowed-sinc filter
%
if (enable_gen_lanczos_winsinc_filter == 1)
    % Setup
    keep enable_gen_lanczos_winsinc_filter enable_gen_least_sq_approx_filter;
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
if (enable_gen_least_sq_approx_filter == 1)
    % Setup
    keep enable_gen_lanczos_winsinc_filter enable_gen_least_sq_approx_filter;
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
