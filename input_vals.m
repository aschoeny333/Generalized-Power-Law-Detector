 % Input Parameters to Run GPL.m as a Script

% Author: Alex Schoeny
% 
% Goal: Provide a Matlab script to easily assign or change values to input
% parameters of GPL.m, run the program, and save the results. See
% descriptions of programs below for meaning of each input parameter. A
% commented line for re-plotting data is also included

% Define input parameters
fnam = 'AU-CZA03-111006-130000.wav';
wav_dir = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/test_wavs';
programs_dir = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev';
pass_band = [150 1800];
stop_band = [100 1850];
t_bounds = [375 450];
gamma = 1;
v1 = 1;
v2 = 2;
eta_thresh = 2.62 * 10^-4 * (75 / (t_bounds(2) - t_bounds(1)))^2;
eta_noise = 2.07 * 10^-5 * (75 / (t_bounds(2) - t_bounds(1)))^2;
t_min = 0.35;
noise_thresh = eta_thresh / 2;
max_noise_dur = 5;
min_noise_dur = 1;
filter_order = 10;
detector_type = 1;
num_receivers = 5;
plot_GPL_data = 0;

% Run GPL.m and save results
[sound, filters, original, whitener_rets, matrices, X_s, intervals, X_masked, ...
    freq_intervals, noise_intervals] = GPL(fnam, wav_dir, programs_dir, pass_band, ...
    stop_band, t_bounds, gamma, v1, v2, eta_thresh, eta_noise, t_min, ...
    noise_thresh, max_noise_dur, min_noise_dur, filter_order, detector_type, plot_GPL_data);

%  Copy the below if re-plotting is necessary, namely with some changes to
%  the variables saved from running GPL.m
% 
%  Plot_Data(original, whitener_rets.m, matrices.N, t_bounds, pass_band, stop_band, ...
%         gamma, v1, v2, eta_thresh, eta_noise, intervals.t, X_masked, ...
%         freq_intervals, noise_intervals.t);

% Define input parameters for associator
rec_dict_tseries = zeros(num_receivers, length(fnam));
rec_dict_tseries(1, :) = fnam;
for i = 2:num_receivers
    next_fnam = fnam;
    next_fnam(8) = num2str(i);
    rec_dict_tseries(i, :) = next_fnam;
end

corr_type = 1;
sig_intervals = intervals.t;
rec_dict_tseries = char(rec_dict_tseries);

% Run associator
[corr_times, range_starts, range_ends] = associator(rec_dict_tseries, corr_type, ...
    sig_intervals, filter_order, freq_intervals, wav_dir, programs_dir);

% Plot associations
Plot_Associations_2(rec_dict_tseries, wav_dir, programs_dir, t_bounds, pass_band, ...
    stop_band, corr_times, sig_intervals, freq_intervals, range_starts, range_ends);