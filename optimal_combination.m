function [optimal_combination_num, optimal_combination_time, optimal_combination_range] = ...
    optimal_combination(remaining_sigs, max_gap, criteria, remaining_freq_intervals, ...
    filter_order, samp_rate, ref_tseries, j_tseries, rec_dict_tseries, rec_j, wav_dir)


    % Determine number of signals to test combining by checking how many
    % adjacent signals are within max_gap seconds of eachother
    num_sigs_to_check = length(remaining_sigs(1,:)) - 1;
    for i = 1:length(remaining_sigs(1,:)) - 1
        if remaining_sigs(1, i+1) - remaining_sigs(2, i) > max_gap
            num_sigs_to_check = i - 1;
            break;
        end
    end

    if num_sigs_to_check > 0
        % Create arrays to save values for each combination, one of which
        % will be the optimal returned value
        criteria_vals = zeros(1, num_sigs_to_check + 1);
        criteria_val_times = zeros(1, num_sigs_to_check + 1);
        range_start_times = zeros(1, num_sigs_to_check + 1);
        range_end_times = zeros(1, num_sigs_to_check + 1);


        % Normalized Cross-Correlation Case
        if criteria == 1
            for i = 0:num_sigs_to_check
                % Create filter for both time series
                cur_freq_intervals = [min(remaining_freq_intervals(1,1), ...
                    remaining_freq_intervals(1, i+1)), max(remaining_freq_intervals(2,1), ...
                    remaining_freq_intervals(2, i+1))];
                bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                    'HalfPowerFrequency1',cur_freq_intervals(1),'HalfPowerFrequency2', ...
                    cur_freq_intervals(2), 'SampleRate', samp_rate);

                % Filter and trim both time series
                ref_tseries_filt = filter(bpFilt, ref_tseries);
                j_tseries_filt = filter(bpFilt, j_tseries);

                % Determine association range for current combination case
                j_duration = possible_range(rec_dict_tseries(1, :), rec_j, [remaining_sigs(1,1), ...
                    remaining_sigs(2, 1 + i)], wav_dir);

                if j_duration(1) < 0
                    j_duration(2) = 0;
                    disp("WARNING: Association range on receiver " + int2str(rec_j) + ", signal " + int2str(i) + " extends before start of audio file, signal detection may fail");
                end
                if j_duration(2) > length(j_tseries_filt) / samp_rate
                    j_duration(2) = length(j_tseries_filt) / samp_rate;
                    disp("WARNING: Association range on receiver " + int2str(rec_j) + ", signal " + int2str(i) + " extends after end of audio file, signal detection may fail");
                end

                range_start_times(1, i+1) = j_duration(1);
                range_end_times(1, i+1) = j_duration(2);

                % Trim time series
                ref_tseries_filt_trimmed = ref_tseries_filt(1 + round(remaining_sigs(1, 1) ...
                    * samp_rate) : round(remaining_sigs(2,1 + i) * samp_rate));
                j_tseries_filt_trimmed = j_tseries_filt(1 + round(j_duration(1) ...
                                * samp_rate) : round(j_duration(2) * samp_rate));

                % Ensure time series are the same length for normalized correlation
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

                % Compute cross correlation
                [r, lags] = xcorr(j_tseries_filt_trimmed, ref_tseries_filt_trimmed, 'normalized');

                % Compute signal to noise ration, approximated by the ratio of
                % maximum to median correlation values
                [max_r, t_lag_1D_ind] = max(r);
                t_lag_1D = lags(t_lag_1D_ind);    

                criteria_vals(i + 1) = max_r;
                criteria_val_times(i + 1) = t_lag_1D;
            end

        % SNR Case
        else
            % Determine the number of 0's to append to each time series so
            % all are the same length as the maximum combination
            % hypothesis case
            j_duration_max = possible_range(rec_dict_tseries(1, :), rec_j, [remaining_sigs(1,1), ...
                remaining_sigs(2, 1 + num_sigs_to_check)], wav_dir);

            if j_duration_max(1) < 0
                j_duration_max(2) = 0;
                disp("WARNING: Association range on receiver " + int2str(rec_j) + ", max combined signal extends before start of audio file, signal detection may fail");
            end
            if j_duration_max(2) > length(j_tseries) / samp_rate
                j_duration_max(2) = length(j_tseries) / samp_rate;
                disp("WARNING: Association range on receiver " + int2str(rec_j) + ", max combined signal extends after end of audio file, signal detection may fail");
            end
            
            max_j_length = (round(j_duration_max(2) * samp_rate)) - ...
                (round(j_duration_max(1) * samp_rate));
            
            max_ref_length = (round(remaining_sigs(2,1 + num_sigs_to_check) * ...
                samp_rate)) - (round(remaining_sigs(1,1) * samp_rate));
            for i = 0:num_sigs_to_check
                % Create filter for both time series
                cur_freq_intervals = [min(remaining_freq_intervals(1,1), ...
                    remaining_freq_intervals(1, i+1)), max(remaining_freq_intervals(2,1), ...
                    remaining_freq_intervals(2, i+1))];
                bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                    'HalfPowerFrequency1',cur_freq_intervals(1),'HalfPowerFrequency2', ...
                    cur_freq_intervals(2), 'SampleRate', samp_rate);

                % Filter and trim both time series
                ref_tseries_filt = filter(bpFilt, ref_tseries);
                j_tseries_filt = filter(bpFilt, j_tseries);

                % Determine association range for current combination case
                j_duration = possible_range(rec_dict_tseries(1, :), rec_j, [remaining_sigs(1,1), ...
                    remaining_sigs(2, 1 + i)], wav_dir);

                if j_duration(1) < 0
                    j_duration(2) = 0;
                    disp("WARNING: Association range on receiver " + int2str(rec_j) + ", signal " + int2str(i) + " extends before start of audio file, signal detection may fail");
                end
                if j_duration(2) > length(j_tseries_filt) / samp_rate
                    j_duration(2) = length(j_tseries_filt) / samp_rate;
                    disp("WARNING: Association range on receiver " + int2str(rec_j) + ", signal " + int2str(i) + " extends after end of audio file, signal detection may fail");
                end

                range_start_times(1, i+1) = j_duration(1);
                range_end_times(1, i+1) = j_duration(2);

                % Trim time series
                ref_tseries_filt_trimmed = ref_tseries_filt(1 + round(remaining_sigs(1, 1) ...
                    * samp_rate) : round(remaining_sigs(2,1 + i) * samp_rate));
                j_tseries_filt_trimmed = j_tseries_filt(1 + round(j_duration(1) ...
                                * samp_rate) : round(j_duration(2) * samp_rate));

                % Ensure time series are the same length as maximum case
                j_length = length(j_tseries_filt_trimmed);
                ref_length = length(ref_tseries_filt_trimmed);

                if j_length < max_j_length
                    append_zeros = zeros(max_j_length - j_length, 1);
                    j_tseries_filt_trimmed = [j_tseries_filt_trimmed; append_zeros];
                end
                if ref_length < max_ref_length
                    append_zeros = zeros(max_ref_length - ref_length, 1);
                    ref_tseries_filt_trimmed = [ref_tseries_filt_trimmed; append_zeros];
                end

                % Compute cross correlati
                [r, lags] = xcorr(j_tseries_filt_trimmed, ref_tseries_filt_trimmed);
%                 figure; hold on;
%                 plot(j_tseries_filt_trimmed, 'g');
%                 plot(ref_tseries_filt_trimmed, 'b');
%                 pause;
%                 [r, lags] = xcorr(j_tseries_filt_trimmed, ref_tseries_filt_trimmed);
%                 if rec_j == 2 & num_sigs_to_check == 2
%                     figure;
%                     plot(lags / samp_rate, r);
%                     title("Correlation for Receiver 2 and " + num2str(i) + " extra detections");
%                 end

                % Find max correlation value
                [max_r, t_lag_1D_ind] = max(r);
                t_lag_1D = lags(t_lag_1D_ind);
                med_r = median(abs(r(r ~= 0)));

                % criteria_vals(i + 1) = max_r;
                criteria_vals(i + 1) = max_r / med_r;
                disp(med_r);
                criteria_val_times(i + 1) = t_lag_1D;
            end
        end

        disp("Criteria Variables for Receiver " + num2str(rec_j));
        disp(criteria_vals);
        % disp(criteria_val_times);

        % Optimal values are all corresponding to maximum value in
        % criteria_vals (highest normalized correlation or highest SNR)
        [~, optimal_combination_num] = max(criteria_vals);
        optimal_combination_time = range_start_times(1,optimal_combination_num) + ...
            criteria_val_times(optimal_combination_num) / samp_rate;
        optimal_combination_range = [range_start_times(optimal_combination_num), ...
            range_end_times(optimal_combination_num)];
        optimal_combination_num = optimal_combination_num - 1;
    else
        % If the next signal is more than max_gap seconds away, return all 0's
        optimal_combination_num = 0;
        optimal_combination_time = 0;
        optimal_combination_range = [0, 0];
    end