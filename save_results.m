% Save Results of GPL Detector and Associator
% 
% Author: Alex Schoeny, Dr. John Spiesberger
% 
% Goal: Provide a Matlab script to easily assign or change values of input
% paramaters, run the program over multiple .wav files in a directory, and
% save the results. See descriptions of programs below for meaning of each
% input parameter.

% Note: wav_dir must include trailing / so they can be used as
% concatenation prefix for reading / loading various files
wav_dir = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/test_wavs/';
programs_dir = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev';
addpath(programs_dir);

% Assign values to detector input parameters
pass_band = [150 1800];
stop_band = [100 1850];
t_bounds = [400 450];
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
associator_type = 2;
num_receivers = 5;
plot_GPL_data = 0;
old_ts = 0; % Keep as false
corr_type = 1;

test_files = dir(fullfile(wav_dir, '*.wav'));
if isempty(test_files)
    disp('No files in wav_dir');
    quit
end

seen_files = zeros(length(test_files), length(wav_dir) + length(test_files(1).name));
num_seen_files = 0;

% Iterate through .wav files in wav_dir
for file = test_files'
    % Assign values to associator input parameters
    rec_dict_tseries = zeros(num_receivers, length(wav_dir) + length(file.name));
    rec_dict_tseries(1, :) = strcat(wav_dir, file.name);
    for i = 2:num_receivers
        next_fnam = file.name;
        next_fnam(8) = num2str(i);
        rec_dict_tseries(i, :) = strcat(wav_dir, next_fnam);
    end
    rec_dict_tseries = char(rec_dict_tseries);
    
    new_file = 1;
    for i = 1:length(seen_files(:, 1))
        for j = 1:num_receivers
            if seen_files(i, :) == rec_dict_tseries(j, :)
                new_file = 0;
            end
        end
    end
    
    if new_file
        disp('**********************');
        disp('Starting New .wav File');
        disp('**********************');
        
        disp('------ Detector ------');
        % Run GPL Detector
        [sound, filters, original, whitener_rets, matrices, X_s, intervals, ...
            X_masked, freq_intervals, noise_intervals] = GPL(file.name, wav_dir, pass_band, ...
            stop_band, t_bounds, gamma, v1, v2, eta_thresh, eta_noise, t_min, ...
            noise_thresh, max_noise_dur, min_noise_dur, filter_order, detector_type, ...
            plot_GPL_data, old_ts);

        disp('----- Associator -----');
        % Run Associator
        [corr_times, range_starts, range_ends] = associator(rec_dict_tseries, ...
            corr_type, intervals.t, filter_order, freq_intervals, wav_dir);

        disp('----- Saving Results -----');
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
        Range_Starts_R1 = range_starts(1,:).';
        Range_Starts_R2 = range_starts(2,:).';
        Range_Starts_R3 = range_starts(3,:).';
        Range_Starts_R4 = range_starts(4,:).';
        Range_Starts_R5 = range_starts(5,:).';
        Range_Ends_R1 = range_ends(1,:).';
        Range_Ends_R2 = range_ends(2,:).';
        Range_Ends_R3 = range_ends(3,:).';
        Range_Ends_R4 = range_ends(4,:).';
        Range_Ends_R5 = range_ends(5,:).';
        Freq_1 = freq_intervals(1,:).';
        Freq_2 = freq_intervals(2,:).';
        Noise_1 = noise_intervals.t(1,:).';
        Noise_2 = noise_intervals.t(2,:).';

        csv_name = "results_" + file.name(1:end - 4) + ".csv";  % - 4 to remove .wav

        results_table = table(Detections_R1, Detections_R1_end, Associations_R2, ...
            Associations_R3, Associations_R4, Associations_R5, Range_Starts_R1, ...
            Range_Starts_R2, Range_Starts_R3, Range_Starts_R4, Range_Starts_R5, ...
            Range_Ends_R1, Range_Ends_R2, Range_Ends_R3, Range_Ends_R4, Range_Ends_R5, ...
            Freq_1, Freq_2, Noise_1, Noise_2, Pass_band_1, Pass_band_2, Stop_band_1, ...
            Stop_band_2, T_bounds_1, T_bounds_2, Gamma, V1, V2, Eta_thresh, Eta_noise, ...
            T_min, Noise_thresh, Max_noise_dur, Min_noise_dur);

        % Save table as .csv file
        writetable(results_table, csv_name);
        
        % Export data in correct format
        save_data_format(corr_times, noise_intervals.t, rec_dict_tseries, ...
            t_bounds, pass_band, stop_band, filter_order, gamma, v1, v2, ...
            noise_thresh, max_noise_dur, min_noise_dur, intervals.t, freq_intervals, ...
            wav_dir, old_ts);

        disp('----- Plotting Associations -----');
        % Plot associations
        Plot_Associations_2(rec_dict_tseries, t_bounds, pass_band, ...
            stop_band, corr_times, intervals.t, freq_intervals, range_starts, range_ends);
%         
%         save_detection_jpgs(rec_dict_tseries, pass_band, stop_band, ...
%             corr_times, intervals.t, freq_intervals, range_starts, range_ends);

        
        % Update looping variables
        seen_files(num_seen_files + 1 : num_seen_files + num_receivers, :) = rec_dict_tseries;
        num_seen_files = num_seen_files + num_receivers;
        
    end
    
    if num_seen_files == length(seen_files(:,1))
        disp('All files seen');
        break
    end
end
