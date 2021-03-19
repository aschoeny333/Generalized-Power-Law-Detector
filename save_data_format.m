function [] = save_data_format(corr_times, noise_intervals, rec_dict_tseries, ...
    t_bounds, pass_band, stop_band, filter_order, gamma, v1, v2, ...
    noise_thresh, max_noise_dur, min_noise_dur, sig_intervals, freq_intervals, wav_dir)

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
    cur_test_stat = test_stat(X, mu, gamma, v1, v2);
    
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
    nrec = num2str(nrec);
    
    fs_mult = 2;
    fs_mult = num2str(fs_mult);
    
    pass_band_exp = zeros(num_receivers, 2);
    for j=1:num_receivers
        pass_band_exp(j, 1) = freq_intervals(2, i);
        pass_band_exp(j, 2) = freq_intervals(1, i);    
    end
    pass_band_exp = num2str(pass_band_exp);
    
    max_dur_emit_call = 140; % Max detection duration on spreadsheet was ~137 - should this be dependent on sound_id?
    max_dur_emit_call = num2str(max_dur_emit_call);
    
    recn = 1:num_receivers;
    recn = num2str(recn);
    
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
    start_of_noise = num2str(start_of_noise);
    
    start_of_signal = start_of_noise;
    start_of_signal(:, 5) = ones(num_receivers, 1) .* (minute + add_min_signal_start);
    start_of_signal(:, 6) = ones(num_receivers, 1) .* second_signal_start;
    start_of_signal = num2str(start_of_signal);
    
    end_of_signal = start_of_noise;
    end_of_signal(:, 5) = ones(num_receivers, 1) .* (minute + add_min_signal_end);
    end_of_signal(:, 6) = ones(num_receivers, 1) .* second_signal_end;
    end_of_signal = num2str(end_of_signal);
    
    end_of_noise = start_of_noise;
    end_of_noise(:, 5) = ones(num_receivers, 1) .* (minute + add_min_noise_end);
    end_of_noise(:, 6) = ones(num_receivers, 1) .* second_noise_end;
    end_of_noise = num2str(end_of_noise);
    
    export = repmat(' ', [6 + 5*num_receivers, 36]);
    export(1) = array_letter;
    export(2, 1:length(sound_id)) = sound_id;
    export(3) = nrec;
    export(4) = fs_mult;
    export(5:(4 + num_receivers), 1:length(pass_band_exp(1,:))) = pass_band_exp;
    export(5+num_receivers, 1:length(max_dur_emit_call)) = max_dur_emit_call;
    export(6+num_receivers, 1:length(recn)) = recn;
    export(7+num_receivers:(6+2*num_receivers), 1:length(start_of_noise(1,:))) = start_of_noise;
    export(7+2*num_receivers:(6+3*num_receivers), 1:length(start_of_signal(1,:))) = start_of_signal;
    export(7+3*num_receivers:(6+4*num_receivers), 1:length(end_of_signal(1,:))) = end_of_signal;
    export(7+4*num_receivers:(6+5*num_receivers), 1:length(end_of_noise(1,:))) = end_of_noise;
    
    export_file_name = convertCharsToStrings(rec_dict_tseries(1,1:end-4));
    export_file_name = export_file_name + "__signal__" + num2str(i) + ".txt";
    writematrix(export, export_file_name);
end



