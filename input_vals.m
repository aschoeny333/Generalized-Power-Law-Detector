fnam = "test_loud_clean_single.wav";
pass_band = [150 1800];
stop_band = [100 1850];
t_bounds = [360 450];
gamma = 1;
v1 = 1;
v2 = 2;
thresh = 10^-8;
eta_thresh = 2.62 * 10^-4;
eta_noise = 2.07 * 10^-5;
t_min = 0.35;

[sound, filters, original, whitener_rets, matrices, X_s, intervals, X_masked, freq_intervals] = ... 
    GPL(fnam, pass_band, stop_band, t_bounds, gamma, v1, v2, thresh, ...
    eta_thresh, eta_noise, t_min);



 Plot_Data(original, whitener_rets.m, matrices.N, matrices.Nt, t_bounds, pass_band, stop_band, ...
        gamma, v1, v2, eta_thresh, eta_noise, intervals.t, X_masked, freq_intervals);
