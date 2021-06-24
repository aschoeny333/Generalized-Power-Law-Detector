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
    
function [corr_times, range_starts, range_ends] = associator(rec_dict_tseries, ...
    corr_type, sig_intervals, filter_order, freq_intervals, wav_dir)

    % Before executing the program, ensure that there are detections to
    % associate
    if ~isempty(sig_intervals(1,:))
        corr_times = zeros(size(rec_dict_tseries, 1), size(sig_intervals, 2));
        range_starts = zeros(size(rec_dict_tseries, 1), size(sig_intervals, 2));
        range_ends = zeros(size(rec_dict_tseries, 1), size(sig_intervals, 2));
        corr_times(1, :) = sig_intervals(1, :);
        range_starts(1, :) = sig_intervals(1, :);
        range_ends(1, :) = sig_intervals(2, :);
        
        % Iterate through the detections on the reference receiver
        for i = 1:length(sig_intervals(1, :))
            disp('Associating a detection');
            % Step 1: Determine interval of investigation on reference receiver
            ref_duration = sig_intervals(:, i);

            % Step 2: Iterate through other (non-reference) receivers 
            for j = 2 :length(rec_dict_tseries(:,1))   % Start at 2 bc 1 is always reference   
                % Step 2.1: Determine interval of investigation on receiver j
                j_duration = possible_range(rec_dict_tseries(1, :), j, ref_duration, wav_dir);
                range_starts(j,i) = j_duration(1);
                range_ends(j,i) = j_duration(2);
                
                % Step 2.2: Determine correlation function
                if corr_type == 1 % 1D, time series correlation
                    % Step 2.2.1: Read in audio time series from the two
                    % receivers
                    [ref_tseries, samp_rate] = audioread(rec_dict_tseries(1,:)); 
                    j_tseries = audioread(rec_dict_tseries(j,:));

                    % Step 2.2.2: Design bandpass filter according to
                    % detected signal frequency bounds
                    cur_freq_intervals = freq_intervals(:, i);
                    bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                        'HalfPowerFrequency1',cur_freq_intervals(1),'HalfPowerFrequency2', ...
                        cur_freq_intervals(2), 'SampleRate', samp_rate);

                    % Step 2.2.3: Filter and trim the reference signal
                    ref_tseries_filt = filter(bpFilt, ref_tseries);
                    ref_tseries_filt = ref_tseries_filt(1 + round(ref_duration(1) ...
                        * samp_rate) : round(ref_duration(2) * samp_rate));

                    % Step 2.2.4: Filter and trim the investigation
                    % interval on receiver j
                    j_tseries_filt = filter(bpFilt, j_tseries);
                    
                    if j_duration(1) < 0
                        j_duration(1) = 0;
                        disp("WARNING: Association range on receiver " + int2str(j) + ", signal " + int2str(i) + " extends before start of audio file, signal detection may fail");
                    end
                    if j_duration(2) > length(j_tseries_filt) / samp_rate
                        j_duration(2) = length(j_tseries_filt) / samp_rate;
                        disp("WARNING: Association range on receiver " + int2str(j) + ", signal " + int2str(i) + " extends after end of audio file, signal detection may fail");
                    end
                   
                    j_tseries_filt = j_tseries_filt(1 + round(j_duration(1) ...
                        * samp_rate) : round(j_duration(2) * samp_rate));

                    % Step 2.2.5: Determine the cross-correlation and
                    % associated lag values
                    j_length = length(j_tseries_filt);
                    ref_length = length(ref_tseries_filt);
                    
                    if j_length > ref_length
                        append_zeros = zeros(j_length - ref_length, 1);
                        ref_tseries_filt = [ref_tseries_filt; append_zeros];
                    end
                    if ref_length > j_length
                        append_zeros = zeros(ref_length - j_length, 1);
                        j_tseries_filt = [j_tseries_filt; append_zeros];
                    end
   
                    [r, lags] = xcorr(j_tseries_filt, ref_tseries_filt, 'normalized');

                    % Step 2.2.6: Find max of r and associated lag - be
                    % careful about order of arguments used in xcorr
                    [~, t_lag_1D_ind] = max(r);
                    t_lag_1D = lags(t_lag_1D_ind);

                    % Step 2.2.7: Store the resulting start bound
                    t_peak = j_duration(1) + t_lag_1D / samp_rate;
                    
                elseif corr_type == 2 % 2D, matrix correlation
                    % Step 2.2.1: Determine the matrices for correlating
                    rec_dict_matrices = []; % Need to write a subroutine to determine these from rec_dict_tseries
                    
                    % Step 2.2.2: Determine the matrix of correlation
                    % coefficients 
                    corr = normxcorr2(rec_dict_matrices(j, :, j_duration), ...
                        single_sig_mask);
                    
                    % Step 2.2.3: Find the max correlation value
                    [~, t_lag_2D] = find(corr==max(corr(:)));
                    
                    % Step 2.2.4: Store the resulting start bound
                    t_peak = j_duration(1) + t_lag_2D;
                    
                elseif corr_type == 3 % Both 1D and 2D correlation
                    %%% 1D CORRELATION %%%
                    % Step 2.2.1: Read in audio time series from the two
                    % receivers
                    [ref_tseries, samp_rate] = audioread(rec_dict_tseries(1,:)); 
                    j_tseries = audioread(rec_dict_tseries(j,:));

                    % Step 2.2.2: Design bandpass filter according to
                    % detected signal frequency bounds
                    cur_freq_intervals = freq_intervals(:, i);
                    bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                        'HalfPowerFrequency1',cur_freq_intervals(1),'HalfPowerFrequency2', ...
                        cur_freq_intervals(2), 'SampleRate', samp_rate);

                    % Step 2.2.3: Filter and trim the reference signal
                    ref_tseries_filt = filter(bpFilt, ref_tseries);
                    ref_tseries_filt = ref_tseries_filt(1 + round(ref_duration(1) ...
                        * samp_rate) : round(ref_duration(2) * samp_rate));

                    % Step 2.2.4: Filter and trim the investigation
                    % interval on receiver j
                    j_tseries_filt = filter(bpFilt, j_tseries);
                    j_tseries_filt = j_tseries_filt(1 + round(j_duration(1) ...
                        * samp_rate) : round(j_duration(2) * samp_rate));

                    % Step 2.2.5: Determine the cross-correlation and
                    % associated lag values
                    [r, lags] = xcorr(j_tseries_filt, ref_tseries_filt);

                    % Step 2.2.6: Find max of r and associated lag - be
                    % careful about order of arguments used in xcorr
                    [~, t_lag_1D_ind] = max(r);
                    t_lag_1D = lags(t_lag_1D_ind);

                    %%% 2D CORRELATION %%%
                    % Step 2.2.7: Determine the matrices for correlating
                    rec_dict_matrices = []; % Need to write a subroutine to determine these from rec_dict_tseries
                    
                    % Step 2.2.8: Determine the matrix of correlation
                    % coefficients 
                    corr = normxcorr2(rec_dict_matrices(j, :, j_duration), ...
                        single_sig_mask);
                    
                    % Step 2.2.9: Find the max correlation value
                    [~, t_lag_2D] = find(corr==max(corr(:)));

                    %%% DETERMINE PROPER LAG %%%
                    % Note: t_peak currently stores average of lags
                    % determined by 1D and 2D correlations. It may be
                    % necessary to change this calculation depending on the
                    % application / dataset
                    
                    % Step 2.2.10: Store the resulting start bound
                    t_peak = j_duration(1) + (t_lag_1D + t_lag_2D) / 2;
                else
                    disp('Error: improper value of corr_type');
                end
    %         Step 2.3: Store association time in corr_times
                corr_times(j, i) = t_peak;
            end
        end   
    else 
        % Return an empty matrix
        corr_times = zeros(length(rec_dict_tseries), 0);
        range_starts = zeros(size(rec_dict_tseries, 1), size(sig_intervals, 2));
        range_ends = zeros(size(rec_dict_tseries, 1), size(sig_intervals, 2));
    end    
end

    