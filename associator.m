% Associator program Draft

% Inputs
%     rec_dict_matrices - Matrix acting as a 'dictionary' associating
%     receiver number and its masked spectrogram
%
%     rec_dict_tseries - Matrix acting as a 'dictionary' associating
%     receiver number and its original time series (trimmed acording to
%     t_bounds)
%
%     corr_type - defines method of correlation
% 
% Iterate through signals on starting receiver

function [corr_bins] = associator(rec_dict_matrices, rec_dict_tseries, corr_type, ...
    signal_intervals, filter_order, freq_bounds, programs_dir, fnam)

    if ~isempty(signal_intervals(1,:))
        corr_bins = zeros(size(rect_dict), size(signal_intervals, 2));
        corr_bins(1, :) = signal_intervals(1, :);
        
        for i = 1:length(signal_intervals(1, :))
            % Step 1: Determine interval of investigation on reference receiver
            ref_duration = signal_intervals(:, i);

            % Step 2: Iterate through other (non-reference) receivers 
            for j = 2 :length(rec_dict)   % Start at 2 bc 1 is always reference   

                % Step 2.1: Determine interval of investigation on receiver j
                j_duration = possible_range(fnam, j, signal_intervals(:, i), programs_dir);

                % TO DO: Implement interval determination program

                % Step 2.2: Determine correlation function
                if corr_type == 1 % 1D, time series correlation
                    ref_tseries = rec_dict_tseries(1);
                    j_tseries = rec_dict_tseries(j);

                    % Design bandpass filter according to detected signal
                    % frequency bounds
                    cur_freq_bounds = freq_bounds(:, i);
                    bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                        'HalfPowerFrequency1',cur_freq_bounds(1),'HalfPowerFrequency2', ...
                        cur_freq_bounds(2), 'SampleRate', samp_rate);

                    % Filter and trim the reference signal
                    ref_tseries_filt = filter(bpFilt, ref_tseries);
                    ref_tseries_filt = ref_tseries_filt(1 + round(ref_duration(1) ...
                        * samp_rate) : round(ref_duration(2) * samp_rate));

                    % Filter and trim the investigation interval on receiver j
                    j_tseries_filt = filter(bpFilt, j_tseries);
                    j_tseries_filt = j_tseries_filt(1 + round(j_duration(1) ...
                        * samp_rate) : round(j_duration(2) * samp_rate));

                    % Determine the cross-correlation and associated lag values
                    [r, lags] = xcorr(j_tseries_filt, ref_tseries_filt);

                    % Find max of r and associated lag - be careful about order
                    % of arguments used in xcorr
                    [~, t_lag_1D_ind] = max(r);
                    t_lag_1D = lags(t_lag_1D_ind);

                    % Store the resulting bounds
                    t_peak = j_duration(1) + t_lag_1D;
                end

                if corr_type == 2 % 2D, matrix correlation
                    corr = normxcorr2(rec_dict_matrices(j, :, j_duration), ...
                        single_sig_mask);
                    [~, t_lag_2D] = find(corr==max(corr(:)));
                    t_peak = j_duration(1) + t_lag_2D;
                end

                if corr_type == 3 % Both 1D and 2D correlation
                    %%% 1D CORRELATION %%%
                    ref_tseries = rec_dict_tseries(1);
                    j_tseries = rec_dict_tseries(j);

                    % Design bandpass filter according to detected signal
                    % frequency bounds
                    cur_freq_bounds = freq_bounds(:, i);
                    bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                        'HalfPowerFrequency1',cur_freq_bounds(1),'HalfPowerFrequency2', ...
                        cur_freq_bounds(2), 'SampleRate', samp_rate);

                    % Filter and trim the reference signal
                    ref_tseries_filt = filter(bpFilt, ref_tseries);
                    ref_tseries_filt = ref_tseries_filt(1 + round(ref_duration(1) ...
                        * samp_rate) : round(ref_duration(2) * samp_rate));

                    % Filter and trim the investigation interval on receiver j
                    j_tseries_filt = filter(bpFilt, j_tseries);
                    j_tseries_filt = j_tseries_filt(1 + round(j_duration(1) ...
                        * samp_rate) : round(j_duration(2) * samp_rate));

                    % Determine the cross-correlation and associated lag values
                    [r, lags] = xcorr(j_tseries_filt, ref_tseries_filt);

                    % Find max of r and associated lag - be careful about order
                    % of arguments used in xcorr
                    [~, t_lag_1D_ind] = max(r);
                    t_lag_1D = lags(t_lag_1D_ind);

                    %%% 2D CORRELATION %%%
                    corr = normxcorr2(rec_dict_matrices(j, :, j_duration), ...
                        single_sig_mask);
                    [~, t_lag_2D] = find(corr==max(corr(:)));


                    %%% DETERMINE PROPER LAG %%%
                    t_peak = j_duration(1) + (t_lag_1D + t_lag_2D) / 2;
                    % How to do this? Currently average 

                end
    %         Step 2.3: Find max of correlation function and store it
                corr_bins(j, i) = t_peak;
            end

        end   
    else 
        % Return an empty matrix
        corr_bins = [];
    end
end
        
% TO DO: Translate time bins stored in loops into seconds

    