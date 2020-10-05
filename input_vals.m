% Input Parameters to Run GPL.m as a Script

% Author: Alex Schoeny
% 
% Goal: Provide a Matlab script to easily assign or change values to input
% parameters of GPL.m, run the program, and save the results. A commented
% line for re-plotting data is also included

% Definte input parameters
fnam = "test_loud_clean.wav";
pass_band = [150 1800];
stop_band = [100 1850];
t_bounds = [0 90];
gamma = 1;
v1 = 1;
v2 = 2;
thresh = 10^-8;
eta_thresh = 2.62 * 10^-4;
eta_noise = 2.07 * 10^-5;
t_min = 0.35;
noise_thresh = eta_thresh;
max_noise_dur = 5;

% Run GPL.m and save results
[sound, filters, original, whitener_rets, matrices, X_s, intervals, X_masked, freq_intervals] = ... 
    GPL(fnam, pass_band, stop_band, t_bounds, gamma, v1, v2, thresh, ...
    eta_thresh, eta_noise, t_min, noise_thresh, max_noise_dur);


%  Copy the below if re-plotting is necessary, namely with some changes to
%  the variables saved from running GPL.m

 Plot_Data(original, whitener_rets.m, matrices.N, matrices.Nt, t_bounds, pass_band, stop_band, ...
        gamma, v1, v2, eta_thresh, eta_noise, intervals.t, X_masked, ...
        freq_intervals, noise_intervals);
