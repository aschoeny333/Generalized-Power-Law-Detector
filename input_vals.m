% Input Parameters to Run GPL.m as a Script

% Author: Alex Schoeny
% 
% Goal: Provide a Matlab script to easily assign or change values to input
% parameters of GPL.m, run the program, and save the results. A commented
% line for re-plotting data is also included

% Definte input parameters
fnam = 'AU-CZA01-111006-130000.wav';
wav_dir = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/test_wavs';
programs_dir = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev';
pass_band = [150 1800];
stop_band = [100 1850];
t_bounds = [375 450];
gamma = 1;
v1 = 1;
v2 = 2;
eta_thresh = 2.62 * 10^-4;
eta_noise = 2.07 * 10^-5;
t_min = 0.35;
noise_thresh = eta_thresh / 2;
max_noise_dur = 5;
min_noise_dur = 1;

% Run GPL.m and save results
[sound, filters, original, whitener_rets, matrices, X_s, intervals, X_masked, ...
    freq_intervals, noise_intervals] = GPL(fnam, wav_dir, programs_dir, pass_band, stop_band, t_bounds, ...
    gamma, v1, v2, eta_thresh, eta_noise, t_min, noise_thresh, max_noise_dur, min_noise_dur);


%  Copy the below if re-plotting is necessary, namely with some changes to
%  the variables saved from running GPL.m

 Plot_Data(original, whitener_rets.m, matrices.N, t_bounds, pass_band, stop_band, ...
        gamma, v1, v2, eta_thresh, eta_noise, intervals.t, X_masked, ...
        freq_intervals, noise_intervals.t);


