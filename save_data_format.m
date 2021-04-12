function [noise_intervals_start, noise_intervals_end] = save_data_format(corr_times, noise_intervals, rec_dict_tseries, ...
    t_bounds, pass_band, stop_band, filter_order, gamma, v1, v2, ...
    noise_thresh, max_noise_dur, min_noise_dur, sig_intervals, freq_intervals, wav_dir, old_ts)

    length_wav_dir = length(wav_dir);
    noise_intervals_start = zeros(size(corr_times));
    noise_intervals_end = zeros(size(corr_times));
    noise_intervals_start(1, :) = noise_intervals(1,:);
    noise_intervals_end(1, :) = noise_intervals(2,:);

    num_receivers = length(rec_dict_tseries(:,1));

    for i=2:num_receivers
        cur_fnam = rec_dict_tseries(i, :);
        [data, samp_rate] = audioread(cur_fnam);

        data_bounded = data(1 + round(t_bounds(1) * samp_rate) : round(t_bounds(2) * samp_rate));

        lpFilt = [];
        bpFilt = [];
        if(pass_band(1)==0)
            % Low pass filter
            lpFilt = designfilt('lowpassfir', 'pass_bandFrequency', pass_band(2),...
                'StopbandFrequency', stop_band(2), 'pass_bandRipple', 1, ...
                'StopbandAttenuation', 35, 'DesignMethod', 'kaiserwin','SampleRate',samp_rate);
            data_filtered = filter(lpFilt, data_bounded);
        else
            % Bandpass filter
            bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
                'HalfPowerFrequency1',stop_band(1),'HalfPowerFrequency2',stop_band(2), ...
                'SampleRate', samp_rate);
            data_filtered = filter(bpFilt, data_bounded);
        end

        filters.l = lpFilt;
        filters.b = bpFilt;

        nfft_val = 2 ^ nextpow2(0.05 * samp_rate);
        [fourier, freqs] = spectrogram(data_filtered, nfft_val, 0.75 * nfft_val, ...
            nfft_val, samp_rate, 'yaxis');

        hz_per_bin = samp_rate / (2 * (length(freqs) - 1)); % Clarification: range of frequencies is samp_rate / 2
        trim_rows = round(pass_band / hz_per_bin); 

        if trim_rows(1) == 0
            fourier_trimmed = fourier(1 + trim_rows(1):trim_rows(2), :);
        else
            fourier_trimmed = fourier(trim_rows(1):trim_rows(2), :);
        end

        X = abs(fourier_trimmed);
        mu = whitener(X);
        disp('Calculating test statistic to determine noise windows');

        if old_ts
            [~, ~, cur_test_stat] = test_stat(X, mu, gamma, v1, v2);
        else
            [~, ~, cur_test_stat] = new_test_stat(X, mu, gamma, v1, v2);
        end

        [~, cols] = size(X);
        sig_intervals_i = zeros(size(sig_intervals));
        sig_intervals_i(1, :) = corr_times(i, :).';
        sig_intervals_i(2, :) = sig_intervals_i(1, :) + sig_intervals(2, :) - sig_intervals(1, :);
        bin_intervals_i = (sig_intervals_i - t_bounds(1)) / ((t_bounds(2) - t_bounds(1)) / cols);
        bin_intervals_i = round(bin_intervals_i);

        noise_intervals_i= noise_bounds(bin_intervals_i, cur_test_stat, cols, ...
            noise_thresh, t_bounds, max_noise_dur, min_noise_dur);

        noise_intervals_start(i, :) = noise_intervals_i.t(1, :);
        noise_intervals_end(i, :) = noise_intervals_i.t(2, :);
    end

    for i=1:length(sig_intervals(1,:))
        array_letter = rec_dict_tseries(1, end - 20);

        sound_id = 'Automatic Detection'; 

        nrec = num_receivers;

        fs_mult = 2;

        pass_band_exp = zeros(num_receivers, 2);
        for j=1:num_receivers
            pass_band_exp(j, 1) = freq_intervals(1, i);
            pass_band_exp(j, 2) = freq_intervals(2, i);    
        end

        max_dur_emit_call = 140; % Max detection duration on spreadsheet was ~137 - should this be dependent on sound_id?

        recn = 1:num_receivers;

        fnam = rec_dict_tseries(1, length_wav_dir+1:end);
        year = 2000 + str2double(fnam(10:11));
        month = str2double(fnam(12:13));
        day = str2double(fnam(14:15));
        hour = str2double(fnam(17:18));
        minute = str2double(fnam(19:20));

        add_min_noise_start = floor(noise_intervals_start(:, i)/60);
        second_noise_start = mod(noise_intervals_start(:, i), 60);

        add_min_signal_start = floor(corr_times(:, i)/60);
        second_signal_start = mod(corr_times(:, i), 60);

        add_min_signal_end = floor((corr_times(:, i) + sig_intervals(2, i) - ...
            sig_intervals(1,i))/60);
        second_signal_end = mod((corr_times(:, i) + sig_intervals(2, i) - ...
            sig_intervals(1,i)), 60);

        add_min_noise_end = floor(noise_intervals_end(:, i)/60);
        second_noise_end = mod(noise_intervals_end(:, i), 60);

        start_of_noise = zeros(num_receivers, 6);
        start_of_noise(:, 1) = ones(num_receivers, 1) * year;
        start_of_noise(:, 2) = ones(num_receivers, 1) * month;
        start_of_noise(:, 3) = ones(num_receivers, 1) * day;
        start_of_noise(:, 4) = ones(num_receivers, 1) * hour;
        start_of_noise(:, 5) = ones(num_receivers, 1) .* (minute + add_min_noise_start);
        start_of_noise(:, 6) = ones(num_receivers, 1) .* second_noise_start;

        start_of_signal = zeros(num_receivers, 6);
        start_of_signal(:, 1) = ones(num_receivers, 1) * year;
        start_of_signal(:, 2) = ones(num_receivers, 1) * month;
        start_of_signal(:, 3) = ones(num_receivers, 1) * day;
        start_of_signal(:, 4) = ones(num_receivers, 1) * hour;
        start_of_signal(:, 5) = ones(num_receivers, 1) .* (minute + add_min_signal_start);
        start_of_signal(:, 6) = ones(num_receivers, 1) .* second_signal_start;

        end_of_signal = zeros(num_receivers, 6);
        end_of_signal(:, 1) = ones(num_receivers, 1) * year;
        end_of_signal(:, 2) = ones(num_receivers, 1) * month;
        end_of_signal(:, 3) = ones(num_receivers, 1) * day;
        end_of_signal(:, 4) = ones(num_receivers, 1) * hour;
        end_of_signal(:, 5) = ones(num_receivers, 1) .* (minute + add_min_signal_end);
        end_of_signal(:, 6) = ones(num_receivers, 1) .* second_signal_end;

        end_of_noise = zeros(num_receivers, 6);
        end_of_noise(:, 1) = ones(num_receivers, 1) * year;
        end_of_noise(:, 2) = ones(num_receivers, 1) * month;
        end_of_noise(:, 3) = ones(num_receivers, 1) * day;
        end_of_noise(:, 4) = ones(num_receivers, 1) * hour;
        end_of_noise(:, 5) = ones(num_receivers, 1) .* (minute + add_min_noise_end);
        end_of_noise(:, 6) = ones(num_receivers, 1) .* second_noise_end;

        export_file_name = convertCharsToStrings(rec_dict_tseries(1,1:end-4));
        export_file_name = export_file_name + "__signal__" + num2str(i) + ".txt";

        txt_file = fopen(export_file_name, 'w');
        fprintf(txt_file, '%1$c\n', array_letter);
        fprintf(txt_file, '%1$s\n', sound_id);
        fprintf(txt_file, '%1$d\n', nrec);
        fprintf(txt_file, '%1$d\n', fs_mult);
        for j = 1:5
            fprintf(txt_file, '%d %d\n', pass_band_exp(j, :));
        end
        fprintf(txt_file, '%d\n', max_dur_emit_call);
        fprintf(txt_file, '%d %d %d %d %d\n', recn);
        for j = 1:5
            fprintf(txt_file, '%1d %1d %1d %1d %1d %1f \n', start_of_noise(j, :));
        end
        for j = 1:5
            fprintf(txt_file, '%1d %1d %1d %1d %1d %1f \n', start_of_signal(j, :));
        end
        for j = 1:5
            fprintf(txt_file, '%1d %1d %1d %1d %1d %1f \n', end_of_signal(j, :));
        end
        for j = 1:5
            fprintf(txt_file, '%1d %1d %1d %1d %1d %1f \n', end_of_noise(j, :));
        end
        fclose(txt_file);
    end