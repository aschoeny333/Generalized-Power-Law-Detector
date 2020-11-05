% Associator program Draft

% Inputs
%     rec_dict_matrices - Some type of "dictionary" associating receiver
%     number and its masked spectrogram rec_dict_tseries - Some type of
%     "dictionary" associating receiver number and its original time series
%     (trimmed acording to t_bounds)
%     corr_type - defines method of correlation
% 
% Iterate through signals on starting receiver
if ~isempty(signal_intervals(1,:))
    corr_bins = zeros(size(rect_dict), size(signal_intervals, 2));
    corr_bins(1, :) = signal_intervals(1, :);
    
    for i = 1:length(signal_intervals(1, :))
%     Step 1: Determine interval of investigation on reference receiver
        ref_duration = signal_intervals(:, i);
    
%     Step 2: Iterate through other (non-reference) receivers 
        for j = 2 :length(rec_dict)   % Start at 2 bc 1 is always reference   
            
%         Step 2.1: Determine interval of investigation on receiver j
            j_duration;
            
            % TO DO: Implement interval determination program

%         Step 2.2: Determine correlation function
            if corr_type == 1 % 1D, time series correlation
                ref_tseries = rec_dict_tseries(1);
                j_tseries = rec_dict_tseries(j);
                
                % Filter the data according to freq_bounds of the reference
                % signal
                cur_freq_bounds = freq_bounds(:, i);
                filter_order=10;
                bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                    'HalfPowerFrequency1',cur_freq_bounds(1),'HalfPowerFrequency2', ...
                    cur_freq_bounds(2), 'SampleRate', samp_rate);
                ref_tseries_filt = filter(bpFilt, ref_tseries);
                j_tseries_filt = filter(bpFilt, j_tseries);
                
                [r, lags] = xcorr(j_tseries_filt, ref_tseries_filt);
                
                % Find max of r and associated lag - be careful about order
                % of arguments used in xcorr
                [max_corr, max_corr_lag_ind] = max(r);
                max_corr_lag = lags(max_corr_lag_ind);
                
                % Store the resulting bounds
                t_peak = j_duration(1) + max_corr_lag;
            end
            
            if corr_type == 2 % 2D, matrix correlation
                corr = normxcorr2(rec_dict_matrices(j, :, j_duration, ...
                    single_sig_mask));
                [~, t_lag] = find(corr==max(corr(:)));
                t_peak = j_duration(1) + t_peak;
            end
            
            if corr_type == 3 % Both 1D and 2D correlation
                %%% 1D CORRELATION %%%
                ref_tseries = rec_dict_tseries(1);
                j_tseries = rec_dict_tseries(j);
                
                % Filter the data according to freq_bounds of the reference
                % signal
                cur_freq_bounds = freq_bounds(:, i);
                filter_order=10;
                bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                    'HalfPowerFrequency1',cur_freq_bounds(1),'HalfPowerFrequency2', ...
                    cur_freq_bounds(2), 'SampleRate', samp_rate);
                ref_tseries_filt = filter(bpFilt, ref_tseries);
                j_tseries_filt = filter(bpFilt, j_tseries);
                
                [r, lags] = xcorr(j_tseries_filt, ref_tseries_filt);
                
                % Find max of r and associated lag - be careful about order
                % of arguments used in xcorr
                [max_corr, max_corr_lag_ind] = max(r);
                max_corr_lag = lags(max_corr_lag_ind);
                
                %%% 2D CORRELATION %%%
                corr = normxcorr2(rec_dict_matrices(j, :, j_duration, ...
                    single_sig_mask));
                [~, t_lag] = find(corr==max(corr(:)));
                
                
                %%% DETERMINE PROPER LAG %%%
                t_peak;
                % How to do this - average? 
                
            end
%         Step 2.3: Find max of correlation function and store it
            corr_bins(j, i) = t_peak;
        end
  
    end   
end
        
% TO DO: Translate time bins stored in loops into seconds

    