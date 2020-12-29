% Save Results of GPL Detector and Associator
% 
% Author: Alex Schoeny
% 
% Goal: Provide a Matlab script to easily assign or change values of input
% paramaters, run the program over multiple .wav files in a directory, and
% save the results. See descriptions of programs below for meaning of each
% input parameter.

wav_dir = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/test_wavs';
programs_dir = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev';

test_files = dir(fullfile(wav_dir, '*.wav'));

% Iterate through .wav files in wav_dir
for file = test_files'
    % Assign values to input parameters
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
    
    % Run GPL Detector
    [sound, filters, original, whitener_rets, matrices, X_s, intervals, ...
        X_masked, freq_intervals, noise_intervals] = GPL(file.name, wav_dir, programs_dir, pass_band, ...
        stop_band, t_bounds, gamma, v1, v2, eta_thresh, eta_noise, t_min, ...
        noise_thresh, max_noise_dur, min_noise_dur, filter_order, detector_type);
    
    % Run associator
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
    
    corr_times = associator(rec_dict_tseries, corr_type, sig_intervals, ...
        filter_order, freq_intervals, wav_dir, programs_dir);
    
    % Save detector outputs and input parameters in a table
    detections = intervals.t;
    num_rows = length(detections(1, :));
    Pass_band_1 = ones(num_rows, 1) * pass_band(1);
    Pass_band_2 = ones(num_rows, 1) * pass_band(2);
    Stop_band_1 = ones(num_rows, 1) * stop_band(1);
    Stop_band_2 = ones(num_rows, 1) * stop_band(2);
    T_bounds_1 = ones(num_rows, 1) * t_bounds(1);
    T_bounds_2 = ones(num_rows, 1) * t_bounds(2);
    Gamma = ones(num_rows, 1) * gamma;
    V1 = ones(num_rows, 1) * v1;
    V2 = ones(num_rows, 1) * v2;
    Eta_thresh = ones(num_rows, 1) * eta_thresh;
    Eta_noise = ones(num_rows, 1) * eta_noise;
    T_min = ones(num_rows, 1) * t_min;
    Noise_thresh = ones(num_rows, 1) * noise_thresh;
    Max_noise_dur = ones(num_rows, 1) * max_noise_dur;
    Min_noise_dur = ones(num_rows, 1) * min_noise_dur;
    Detections_R1 = detections(1,:).';
    Detections_R1_end = detections(2,:).';
    Associations_R2 = corr_times(2, :).';
    Associations_R3 = corr_times(3, :).';
    Associations_R4 = corr_times(4, :).';
    Associations_R5 = corr_times(5, :).';
    Freq_1 = freq_intervals(1,:).';
    Freq_2 = freq_intervals(2,:).';
    Noise_1 = noise_intervals.t(1,:).';
    Noise_2 = noise_intervals.t(2,:).';
    
    csv_name = "results_" + file.name(1:end - 4) + ".csv";  % - 4 to remove .wav
    
    results_table = table(Detections_R1, Detections_R1_end, Associations_R2, ...
        Associations_R3, Associations_R4, Associations_R5, Freq_1, Freq_2, Noise_1, ...
        Noise_2, Pass_band_1, Pass_band_2, Stop_band_1, Stop_band_2, T_bounds_1, ...
        T_bounds_2, Gamma, V1, V2, Eta_thresh, Eta_noise, T_min, Noise_thresh, ...
        Max_noise_dur, Min_noise_dur);
    
    % Save table as .csv file
    writetable(results_table, csv_name);
end
