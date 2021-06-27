% Plot Correlations of a Signal Between Two Receivers
% 
% Author: Alex Schoeny, Dr. John Spiesberger
% 
% Goal: For a specific signal s and a specific receiver r, plot the
% correlation values for times on r. Useful for debugging in associator.m.
% 
% Inputs:
% 
%     receiver - 1 x 1 double, integer value for which receiver to
%     correlate with
%     
%     signal - 1 x 1 double, integer value for which detection to correlate
%     with
%   
%     sig_intervals - Double matrix of size 2 x k'' as defined in
%     detector.m, where each column is of the form [start; end], where
%     start and end define the beginning and ending time value in seconds
%     of the signal identified by the detector
%
%     rec_dict_tseries - Char array of size r x l, with row i corresponding
%     to the name of the .wav file for receiver i.
%
%     freq_intervals - Double matrix of size 2 x k'' as defined in
%     detector.m, where each column is of the form [start; end], where
%     start and end define the beginning and ending frequency bins of the
%     signal contained in X_masked. Useful as input to associator
%
%     fnam - Char array, .wav file name (e.g. 'test.wav')
% 
%     programs_dir - Char array, gives directory containing associator.m
%     and other relevant programs
%
%     wav_dir - Char array, gives directory containing .wav files
%
%     filter_order - 1 x 1 Double, input argument for Matlab function
%     designfilt, suggested value of 10
%
% Outputs:
% 
%     test_j_duration - 1 x 2 Double, time bounds on non-reference receiver
%     which the detection in question from reference receiver must occur.
%     Calculated by possible_range.m

function [test_j_duration] = Plot_Correlations(receiver, signal, sig_intervals, ...
    rec_dict_tseries, freq_intervals, fnam, programs_dir, wav_dir, filter_order)

    test_ref_duration = sig_intervals(:, signal);
    test_j_duration = possible_range(fnam, receiver, test_ref_duration, programs_dir);
    test_j_duration(1) = test_j_duration(1) - (test_ref_duration(2) - test_ref_duration(1));
    % test_j_duration = [430, 450];

    cd(wav_dir)
    [ref_tseries, samp_rate] = audioread(rec_dict_tseries(1, :)); % : is important - elts are char arrays
    j_tseries = audioread(rec_dict_tseries(5,:));
    cd(programs_dir)

    % Design bandpass filter according to detected signal
    % frequency bounds
    cur_freq_intervals = freq_intervals(:, 1);
    bpFilt = designfilt('bandpassiir','FilterOrder',filter_order, ...
        'HalfPowerFrequency1',cur_freq_intervals(1),'HalfPowerFrequency2', ...
        cur_freq_intervals(2), 'SampleRate', samp_rate);

    % Filter and trim the reference signal
    ref_tseries_filt = filter(bpFilt, ref_tseries);
    ref_tseries_filt = ref_tseries_filt(1 + round(test_ref_duration(1) ...
        * samp_rate) : round(test_ref_duration(2) * samp_rate));

    % Filter and trim the investigatsion interval on receiver j
    j_tseries_filt = filter(bpFilt, j_tseries);
    j_tseries_filt = j_tseries_filt(1 + round(test_j_duration(1) ...
        * samp_rate) : round(test_j_duration(2) * samp_rate));

    % Determine the cross-correlation and associated lag values
    ref_tseries_filt_longer = zeros(size(j_tseries_filt));
    ref_tseries_filt_longer(1:length(ref_tseries_filt)) = ref_tseries_filt;
    [r, lags] = xcorr(j_tseries_filt, ref_tseries_filt_longer, 'biased');

    figure;
    plot(test_j_duration(1) + lags / samp_rate, r);
    title("Correlations of receiver " + num2str(receiver) + " and signal " + num2str(signal));
    xlabel("Time on receiver " + num2str(receiver));
    ylabel("Correlation value");
    
end