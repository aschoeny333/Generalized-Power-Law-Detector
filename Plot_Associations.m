% Plot Associations of All Detections Across a Set of Receivers
% 
% Author: Alex Schoeny
% 
% Goal: Plot the time bounds of all detections on a reference receiver, as
% well as plotting the associations of each of those detections determined
% by associator.m on a set of other receivers
% 
% Inputs:
%
%     rec_dict_tseries - Char array of size r x l, with row i corresponding
%     to the name of the .wav file for receiver i.
%
%     wav_dir - Char array, gives directory containing .wav files
% 
%     programs_dir - Char array, gives directory containing associator.m
%     and other relevant programs
%
%     t_bounds - 1 x 2 Double array, time interval of interest of form [start
%     time, end time] in sec (e.g. [360, 450]). 
%
%     pass_band - 1 x 2 Double array, Min/max passband frequencies (Hz).
%     Minimum in pass_band(1). For a lowpass filter, set pass_band(1)=0. If
%     pass_band is empty, output data are not filtered.
% 
%     stop_band - 1 x 2 Double array, Min/max stop band (Hz) for bandpass
%     filter. Minimum in stop_band(1). stop_band(1) only used when
%     pass_band(1)>0. then stop_band(1) <pass_band(1). Only used when
%     pass_band is not empty.
%
%     corr_times - Double matrix of size r x k'' as defined as detector.m,
%     where entry i,j corresponds to the start time of signal j on receiver
%     i. End times can be determined by adding corresponding duration of
%     the detection from the reference receiver, which can be calculated
%     from sig_intervals
%
%     sig_intervals - Double matrix of size 2 x k'' as defined in
%     detector.m, where each column is of the form [start; end], where
%     start and end define the beginning and ending time value in seconds
%     of the signal identified by the detector

function [] = Plot_Associations(rec_dict_tseries, wav_dir, programs_dir, ...
    t_bounds, pass_band, stop_band, corr_times, sig_intervals, freq_intervals, ...
    range_starts, range_ends)

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
%             xline(corr_times(i, 1), 'color', 'r', 'LineWidth', 2);
%             xline(corr_times(i, 1) + sig_intervals(2,1) - sig_intervals(1,1), ...
%                 'color', 'r', 'LineWidth', 2);
    %         yline(freq_intervals(1, 1), 'color', 'g', 'LineWidth', 2);
    %         yline(freq_intervals(2, 1), 'color', 'r', 'LineWidth', 2);
%             rectangle('Position', [corr_times(i,1), freq_intervals(1,1), ...
%                 sig_intervals(2,1) - sig_intervals(1,1), freq_intervals(2,1) - freq_intervals(1,1)], ...
%                 'EdgeColor', 'r', 'LineWidth', 2); 
            verts = [corr_times(i,1) freq_intervals(1, 1) 5; corr_times(i,1) freq_intervals(2, 1) 5;...
                corr_times(i,1)+(sig_intervals(2,1)-sig_intervals(1,1)) freq_intervals(2, 1) 5; ...
                corr_times(i,1)+(sig_intervals(2,1)-sig_intervals(1,1)) freq_intervals(1, 1) 5]; 
            face = [1 2 3 4];
            patch('Faces', face, 'Vertices', verts, 'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 2);
            
            if i ~= 1
                xline(range_starts(i,1), '--', 'color', 'r', 'LineWidth', 2);
                xline(range_ends(i,1), '--', 'color', 'r', 'LineWidth', 2);
            end
        else
            for j = 1 : length(corr_times(1, :))
%                 xline(corr_times(i, j), 'color', 'r', 'LineWidth', 2);
%                 xline(corr_times(i, j) + sig_intervals(2,j) - sig_intervals(1,j), ...
%                     'color', 'r', 'LineWidth', 2);
    %             yline(freq_intervals(1, i), 'color', 'g', 'LineWidth', 2);
    %             yline(freq_intervals(2, i), 'color', 'r', 'LineWidth', 2);
%                 rectangle('Position', [corr_times(i, j), freq_intervals(1,j), ...
%                 sig_intervals(2,j) - sig_intervals(1,j), freq_intervals(2,j) - freq_intervals(1,j)], ...
%                 'EdgeColor', 'r', 'LineWidth', 2); 
                verts = [corr_times(i,j) freq_intervals(1, j) 5; corr_times(i,j) freq_intervals(2, j) 5;...
                    corr_times(i,j)+(sig_intervals(2,j)-sig_intervals(1,j)) freq_intervals(2, j) 5; ...
                    corr_times(i,j)+(sig_intervals(2,j)-sig_intervals(1,j)) freq_intervals(1, j) 5]; 
                face = [1 2 3 4];
                color_list = ['r', 'g', 'y', 'w', 'c'];
                patch('Faces', face, 'Vertices', verts, 'EdgeColor', color_list(mod(j,5)), ...
                    'FaceColor', 'none', 'LineWidth', 1);
                
                if i ~= 1
                    xline(range_starts(i,j), '--', 'color', color_list(mod(j, 5)), 'LineWidth', 2);
                    xline(range_ends(i,j), '--', 'color', color_list(mod(j, 5)), 'LineWidth', 2);
                end
            end
        end
        xlabel('Time (secs)');
        c = colorbar;
        ylabel('Frequency(HZ)');
        ylabel(c, 'Power');
        title("Associated signals on receiver " + i);
    end
end