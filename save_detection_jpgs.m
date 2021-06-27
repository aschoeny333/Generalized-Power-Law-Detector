% Save .jpg files of the spectrogram of each detection
% 
% Author: Alex Schoeny, Dr. John Spiesberger
% 
% Goal: Save an image of the spectrogram of each detection across all
% receivers after the execution of the software package across a given set
% of audio files

function [] = save_detection_jpgs(rec_dict_tseries, pass_band, stop_band, ...
    corr_times, sig_intervals, freq_intervals, range_starts, range_ends)

    for i = 1:length(sig_intervals(1,:))
        disp('Saving jpeg of a detection');
        figure;
        earliest_time = min(range_starts(:, i));
        latest_time = max(range_ends(:, i));

        for j = 1:length(rec_dict_tseries(:,1))
            subplot(length(rec_dict_tseries(:,1)), 1, j)
            [data, samp_rate] = audioread(rec_dict_tseries(j, :));
            
            if earliest_time < 1
                earliest_time = 1;
            end
            if latest_time > length(data) / samp_rate - 1
                latest_time = length(data) / samp_rate - 1;
            end

            data_bounded = data(1 + round((earliest_time - 1) * samp_rate) : ...
                round((latest_time + 1) * samp_rate));

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

            line_height = max(max(X));

            surf(times + earliest_time - 1, freqs_trimmed, X, 'EdgeColor', 'none');
            axis xy; 
            axis tight; 
            colormap(jet); view(0,90);

            color_list = ['r', 'g', 'y', 'w', 'c'];
            verts = [corr_times(j,i) freq_intervals(1, i) line_height; corr_times(j,i) freq_intervals(2, i) line_height;...
                corr_times(j,i)+(sig_intervals(2,i)-sig_intervals(1,i)) freq_intervals(2, i) line_height; ...
                corr_times(j,i)+(sig_intervals(2,i)-sig_intervals(1,i)) freq_intervals(1, i) line_height]; 
            face = [1 2 3 4];
            patch('Faces', face, 'Vertices', verts, 'EdgeColor', color_list(mod(i,5) + 1), 'FaceColor', 'none', 'LineWidth', 0.5);

            if j > 1
                xline(range_starts(j, i), '--', 'color', color_list(mod(i, 5) + 1), 'LineWidth', 2);
                xline(range_ends(j, i), '--', 'color', color_list(mod(i, 5) + 1), 'LineWidth', 2);
            end

            xlabel('Time (secs)');
            c = colorbar;
            ylabel('Frequency(HZ)');
            ylabel(c, 'Power');
            title("Associated signals on receiver 1");

        end
        save_name = strcat(rec_dict_tseries(1, end-25:end-4), '_signal_', num2str(i), '.jpg');
        vers_string = version;
        if str2double(vers_string(17:20)) > 2019
            exportgraphics(gcf, save_name, 'Resolution', 500);
        else
            saveas(gcf, save_name);
        end
        clf; close;
    end