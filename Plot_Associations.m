function [] = Plot_Associations(rec_dict_tseries, wav_dir, programs_dir, ...
    t_bounds, pass_band, stop_band, corr_times, signal_intervals)

    figure; 
    num_receivers = length(rec_dict_tseries(:,1));
    
    for i = 1:num_receivers
        cd(wav_dir);
        [data, samp_rate] = audioread(rec_dict_tseries(i, :));
        cd(programs_dir);
        
        data_bounded = data(1 + round(t_bounds(1) * samp_rate) : ...
            round(t_bounds(2) * samp_rate));

        if(pass_band(1)==0)
            % Low pass filter
            lpFilt = designfilt('lowpassfir', 'pass_bandFrequency', pass_band(2),...
                'StopbandFrequency', stop_band(2), 'pass_bandRipple', 1, ...
                'StopbandAttenuation', 35, 'DesignMethod', 'kaiserwin','SampleRate',samp_rate);
            data_filtered = filter(lpFilt, data_bounded);
        else
            % Bandpass filter
            filter_order=10;
            bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                'HalfPowerFrequency1',stop_band(1),'HalfPowerFrequency2',stop_band(2), ...
                'SampleRate', samp_rate);
            data_filtered = filter(bpFilt, data_bounded);
        end

        nfft_val = 2 ^ nextpow2(0.05 * samp_rate);
        [fourier, freqs, times] = spectrogram(data_filtered, nfft_val, 0.75 * nfft_val, ...
            nfft_val, samp_rate, 'yaxis');

        hz_per_bin = samp_rate / (2 * (length(freqs) - 1)); % Clarification: range of frequencies is samp_rate / 2
        trim_rows = round(pass_band / hz_per_bin); 
        if trim_rows(1) == 0
            freqs_trimmed = freqs(1 + trim_rows(1):trim_rows(2));
            fourier_trimmed = fourier(1 + trim_rows(1):trim_rows(2), :);
        else
            freqs_trimmed = freqs(trim_rows(1):trim_rows(2));
            fourier_trimmed = fourier(trim_rows(1):trim_rows(2), :);
        end
        
        X = abs(fourier_trimmed);
        
        subplot(num_receivers, 1, i)
        surf(times + t_bounds(1), freqs_trimmed, X, 'EdgeColor', 'none');
        axis xy; 
        axis tight; 
        colormap(jet); view(0,90);
        if numel(corr_times) == 0
            % Don't plot any time bound lines
        elseif length(corr_times(1, :)) == 1
            xline(corr_times(i, 1), 'color', 'r', 'LineWidth', 2);
            xline(corr_times(i, 1) + signal_intervals(2,1) - signal_intervals(1,1), ...
                'color', 'r', 'LineWidth', 2);
    %         yline(freq_intervals(1, 1), 'color', 'g', 'LineWidth', 2);
    %         yline(freq_intervals(2, 1), 'color', 'r', 'LineWidth', 2);
    %         rectangle('Position', [time_intervals(1,1), freq_intervals(1,1), ...
    %             time_intervals(2,1) - time_intervals(1,1), freq_intervals(2,1) - freq_intervals(1,1)], ...
    %             'EdgeColor', 'w', 'LineWidth', 3); 
        else
            for j = 1 : length(corr_times(1, :))
                xline(corr_times(i, j), 'color', 'r', 'LineWidth', 2);
                xline(corr_times(i, j) + signal_intervals(2,j) - signal_intervals(1,j), ...
                    'color', 'r', 'LineWidth', 2);
    %             yline(freq_intervals(1, i), 'color', 'g', 'LineWidth', 2);
    %             yline(freq_intervals(2, i), 'color', 'r', 'LineWidth', 2);
    %             rectangle('Position', [time_intervals(1,i), freq_intervals(1,i), ...
    %             time_intervals(2,i) - time_intervals(1,i), freq_intervals(2,i) - freq_intervals(1,i)], ...
    %             'EdgeColor', 'k', 'LineWidth', 3); 
            end
        end
        xlabel('Time (secs)');
        c = colorbar;
        ylabel('Frequency(HZ)');
        ylabel(c, 'Power');
        title("Associated signals on receiver " + i);
    end
end