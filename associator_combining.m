% Signal Associator
% 
% Author: Alex Schoeny
% 
% Goal: For each detection on a reference receiver, loop through other
% receivers and use cross correlation to determine the start time of the
% same signal
% 
% Inputs:
%
%     rec_dict_tseries - Char array of size r x l, with row i corresponding
%     to the name of the .wav file for receiver i.
% 
%     corr_type - 1 x 1 Double, defines method of correlation. If 1,
%     correlates the time series; if 2, correlates the fourier matrices; if
%     3, correlates both and takes average of correlated time bins
% 
%     sig_intervals - Double matrix of size 2 x k'' as defined in
%     detector.m, where each column is of the form [start; end], where
%     start and end define the beginning and ending time value in seconds
%     of the signal identified by the detector
% 
%     filter_order - 1 x 1 Double, input argument for Matlab function
%     designfilt, suggested value of 10
%     
%     freq_intervals - Double matrix of size 2 x k'' as defined in
%     detector.m, where each column is of the form [start; end], where
%     start and end define the beginning and ending frequency bins of the
%     signal contained in X_masked. Useful as input to associator
% 
%     wav_dir - Char array, gives directory containing .wav files
% 
%     programs_dir - Char array, gives directory containing associator.m
%     and other relevant programs
%
% Outputs:
% 
%     corr_times - Double matrix of size r x k'' as defined as detector.m,
%     where entry i,j corresponds to the start time of signal j on receiver
%     i. End times can be determined by adding corresponding duration of
%     the detection from the reference receiver, which can be calculated
%     from sig_intervals
    
function [corr_times, range_starts, range_ends, new_sig_intervals, new_freq_intervals, ...
    new_noise_intervals] = associator_combining(rec_dict_tseries, ...
    corr_type, sig_intervals, filter_order, freq_intervals, wav_dir, noise_intervals, ...
    max_gap, criteria, combine_rule)

    % Before executing the program, ensure that there are detections to
    % associate
    if ~isempty(sig_intervals(1,:))
        % Set up return variables
        corr_times = zeros(size(rec_dict_tseries, 1), size(sig_intervals, 2));
        range_starts = zeros(size(rec_dict_tseries, 1), size(sig_intervals, 2));
        range_ends = zeros(size(rec_dict_tseries, 1), size(sig_intervals, 2));
        corr_times(1, :) = sig_intervals(1, :);
        range_starts(1, :) = sig_intervals(1, :);
        range_ends(1, :) = sig_intervals(2, :);
        
        % Iterate through the detections on the reference receiver
        for i = 1:length(sig_intervals(1, :))
            disp('Associating a detection');
            % Determine interval of investigation on reference receiver
            ref_duration = sig_intervals(:, i);
            
            % Set up combining variables - used to determine if group of
            % adjacent signals ought to be combined
            optimal_combination_nums = zeros(1, length(rec_dict_tseries(:, 1)) - 1);
            optimal_combination_times = zeros(1, length(rec_dict_tseries(:, 1)) - 1);
            optimal_combination_starts = zeros(1, length(rec_dict_tseries(:, 1)) - 1);
            optimal_combination_ends = zeros(1, length(rec_dict_tseries(:, 1)) - 1);

            % Iterate through other (non-reference) receivers 
            for j = 2 :length(rec_dict_tseries(:,1))   % Start at 2 bc 1 is always reference   
                % Determine interval of investigation on receiver j
                j_duration = possible_range(rec_dict_tseries(1, :), j, ref_duration, wav_dir);
                range_starts(j,i) = j_duration(1);
                range_ends(j,i) = j_duration(2);
                
                % Determine correlation function
                if corr_type == 1 % 1D, time series correlation
                    % Read in audio time series from the two receivers
                    [ref_tseries, samp_rate] = audioread(rec_dict_tseries(1,:)); 
                    j_tseries = audioread(rec_dict_tseries(j,:));

                    % Design bandpass filter according to detected signal
                    % frequency bounds
                    cur_freq_intervals = freq_intervals(:, i);
                    bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                        'HalfPowerFrequency1',cur_freq_intervals(1),'HalfPowerFrequency2', ...
                        cur_freq_intervals(2), 'SampleRate', samp_rate);

                    % Filter and trim the reference signal
                    ref_tseries_filt = filter(bpFilt, ref_tseries);
                    ref_tseries_filt_trimmed = ref_tseries_filt(1 + round(ref_duration(1) ...
                        * samp_rate) : round(ref_duration(2) * samp_rate));

                    % Filter and trim the investigation interval on
                    % receiver j
                    j_tseries_filt = filter(bpFilt, j_tseries);
                    
                    if j_duration(1) < 0
                        j_duration(2) = 0;
                        disp("WARNING: Association range on receiver " + int2str(j) + ", signal " + int2str(i) + " extends before start of audio file, signal detection may fail");
                    end
                    if j_duration(2) > length(j_tseries_filt) / samp_rate
                        j_duration(2) = length(j_tseries_filt) / samp_rate;
                        disp("WARNING: Association range on receiver " + int2str(j) + ", signal " + int2str(i) + " extends after end of audio file, signal detection may fail");
                    end
                   
                    j_tseries_filt_trimmed = j_tseries_filt(1 + round(j_duration(1) ...
                        * samp_rate) : round(j_duration(2) * samp_rate));

                    % Determine the cross-correlation and associated lag
                    % values after appending zeros if necessary to ensure
                    % vectors are of equal lengths
                    j_length = length(j_tseries_filt_trimmed);
                    ref_length = length(ref_tseries_filt_trimmed);
                    
                    if j_length > ref_length
                        append_zeros = zeros(j_length - ref_length, 1);
                        ref_tseries_filt_trimmed = [ref_tseries_filt_trimmed; append_zeros];
                    end
                    if ref_length > j_length
                        append_zeros = zeros(ref_length - j_length, 1);
                        j_tseries_filt_trimmed = [j_tseries_filt_trimmed; append_zeros];
                    end
   
                    [r, lags] = xcorr(j_tseries_filt_trimmed, ref_tseries_filt_trimmed, 'normalized');

                    % Find max of r and associated lag time
                    [~, t_lag_1D_ind] = max(r);
                    t_lag_1D = lags(t_lag_1D_ind);

                    % Step 2.2.7: Store the resulting start bound
                    t_peak = j_duration(1) + t_lag_1D / samp_rate;    
                end

                % Run optimal_combination on the current signal - returns
                % number of signals to combine and corresponding start time
                % and correlation range for the given signal and receiver
                [optimal_combination_num, optimal_combination_time, optimal_combination_range] = ...
                    optimal_combination(sig_intervals(:, i:end), max_gap, criteria, ...
                    freq_intervals(:, i:end), filter_order, samp_rate, ...
                    ref_tseries, j_tseries, rec_dict_tseries, j, wav_dir);
                
                % Save results from above
                optimal_combination_nums(j - 1) = optimal_combination_num;
                optimal_combination_times(j - 1) = optimal_combination_time;
                optimal_combination_starts(j - 1) = optimal_combination_range(1);
                optimal_combination_ends(j - 1) = optimal_combination_range(2);
                
                % Store association time in corr_times
                corr_times(j, i) = t_peak;
            end
            
            disp("Signals determined to combine on each receiver");
            disp(optimal_combination_nums);
            
            % Determine number of signals to combine on all receivers based
            % on combine_rule, optimal_combination_nums
            if combine_rule == 1
                combination_count = min(optimal_combination_nums);
            elseif combine_rule == 2
                [combination_count, num_cc] = mode(optimal_combination_nums);
                if num_cc < 2
                    combination_count = min(optimal_combination_nums);
                end
            else
                [combination_count, num_cc] = mode(optimal_combination_nums);
                if num_cc < 3
                    combination_count = min(optimal_combination_nums);
                end
            end
            
            if combination_count > 0
                % Write program to combine everything and return adjusted
                % variables
                corr_times(2:end, i) = optimal_combination_times.';
                corr_times(:, i+1:i+combination_count) = [];
                
                range_starts(2:end, i) = optimal_combination_starts.';
                range_starts(:, i+1:i+combination_count) = [];
                
                range_ends(2:end, i) = optimal_combination_ends.';
                range_ends(:, i+1:i+combination_count) = [];
                
                sig_intervals(:, i) = [sig_intervals(1, i); sig_intervals(2, i+combination_count)];
                sig_intervals(:, i+1:i+combination_count) = [];
                
                freq_intervals(:, i) = [min(freq_intervals(1,i), freq_intervals(1, i+combination_count)); ...
                    max(freq_intervals(2,1), freq_intervals(2, i+combination_count))];
                freq_intervals(:, i+1:i+combination_count) = [];
                
                noise_intervals(:, i) = [noise_intervals(1, i); noise_intervals(2, i + combination_count)];
                noise_intervals(:, i+1:i+combination_count) = [];
            end
            
            if i == length(sig_intervals)
                break;
            end     
        end   
    else 
        % Return an empty matrix
        corr_times = zeros(length(rec_dict_tseries), 0);
        range_starts = zeros(size(rec_dict_tseries, 1), size(sig_intervals, 2));
        range_ends = zeros(size(rec_dict_tseries, 1), size(sig_intervals, 2));
    end  
    
    % Save adjusted variables to return
    new_freq_intervals = freq_intervals;
    new_sig_intervals = sig_intervals;
    new_noise_intervals = noise_intervals;
end

    