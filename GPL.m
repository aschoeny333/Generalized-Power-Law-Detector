% Generalized Power Law Algorithm Implementation 
% 
% Author: Alex Schoeny
% 
% Goal: Produce spectrograms with whitened background noise without losing
% relative power and increasing visual clarity on spectrogram of sounds
% from marine biology
% 
% References:
% Helble, Tyler A et al. ?A generalized power-law detection algorithm for
% humpback whale vocalizations.? The Journal of the Acoustical Society of
% America vol. 131,4 (2012): 2682-99. doi:10.1121/1.3685790
% 
% Inputs:
% 
%     fnam - String, .wav file name (e.g. 'test.wav')
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
%     t_bounds - 1 x 2 Double array, time interval of interest of form [start
%     time, end time] in sec (e.g. [360, 450]). If start time > end time,
%     values are swapped. If any other improper entry (e.g. start time =
%     end time, start time < 0, end time > time series length), t_bounds is
%     set to [0, time series length]
% 
%     gamma - 1 x 1 Double, paramater used to determine A and B (e.g. 1) as
%     used in Helble et al. (2012) equations 7 and 8, p. 2685
% 
%     v1 - 1 x 1 Double, exponential parameter applied to A to determine N
%     (e.g. 1) as defined in Helble et al. (2012) equation 6, p. 2685
% 
%     v2 - 1 x 1 Double, exponential parameter applied to B to determine N
%     (e.g. 2)as defined in Helble et al. (2012) equation 6, p. 2685
% 
%     eta_thresh - 1 x 1 Double, minimum test statistic threshold a signal
%     must reach in at least one time bin in order to be identified by
%     detector.m. Derived largely in sections IIIB-D in Helble et al.
%     (2012) Exact value given in Section IV (2.62 * 10^-4), Helble et al.
%     (2012) p. 2690
% 
%     eta_noise - 1 x 1 Double, minimum test statistic threshold a signal
%     must reach in all time bins in order to have its time interval
%     determined by detector.m; that is, test statistics below this value
%     imply a signal does not exist in that time bin. Derived largely in
%     sections IIIB-D in Helble et al. Exact value given in Section IV
%     (2.07 * 10^-5), Helble et al. (2012) p. 2690
%
%     noise_thresh - 1 x 1 double, value of test statistic above which
%     noise interval iteration procedure stops. Not referenced in either
%     Helble paper, but consider using eta_noise
%
%     max_noise_dur - 1 x 1 double, value in seconds that the length of the
%     noise bounds cannot exceed
%  
%     min_noise_dur - 1 x 1 double, value in seconds that the length of the
%     noise bounds must exceed
%
%     filter_order - 1 x 1 double, input argument for Matlab function
%     designfilt, suggested value of 10
%
%     detector_type - 1 x 1 double. If 1, runs detector.m, which replaces
%     each individual detection with columns containing noise conditions
%     and then recalculates test statistic (hence is better suited for data
%     with many signals, has longer runtime). If any other value, runs
%     detector_simple.m, which does not replace columns of detections
%
% Outputs
%
%     A note on size definitions: Some matrices and arrays sizes are
%     defined using m, n, and m'. The values of these will change based on
%     the values of the inputs above, especially pass_band, stop_band, and
%     t_bounds. Note also that m' <= m. Matrices all have rows with
%     constant frequency and columns with constant time.
%
%     sound - Structured array with fields d, s defined as follows:
%         d - Variable name: data_bounded. Double column vector of size q,
%         where q = samp_rate * bounds_len, the number of acoustic samples;
%         raw numeric data output from audioread(fnam) trimmed according to
%         t_bounds
% 
%         s - Variable name: samp_rate. 1 x 1 Double, sample rate output
%         from audioread in Hz
% 
%     filters - Structured array with fields l, p defined as follows:
%         l - Variable name: lpFilt. lowpass filter. Structured array
%         outputted by designfilt.m. Only meaningful when pass_band(1)==0.
%         Output data are made using d_filt=filter(lpFilt,d). When
%         pass_band(1)=0, it is empty.
% 
%         b - Variable name: bpFilt. bandpass filter. Structured array
%         outputted by designfilt.m. Only meaningful when pass_band(1)>0.
%         Output data are made using d_filt=filter(bpFilt,d). When
%         pass_band(1)>0, it is empty.
%         
%     original - Structured array with fields X, f, t, n defined as
%     follows:
%         X - Variable name: X. Double matrix of size m x n as described in
%         the note above; absolute value of fourier_trimmed. Not identical
%         to the X defined in Helble et al. (2012) p. 2684, but the X in
%         Helble et al (2012) is only ever used in a context where its
%         absolute value is taken, hence the shorthand used here
% 
%         f - Variable name: freqs_trimmed. Double column vector of size m
%         as described in note above; subarray of freqs determined by
%         trim_rows. Values in Hz
% 
%         t - Variable name: times. Double row vector of size n as
%         described in note above; times output from spectrogram function
%         defining the time corresponding to each column of fourier. Values
%         in sec
% 
%         n - Variable name: nfft_val. 1 x 1 Double, 0.05 * samp_rate
%         rounded up to nearest power of 2, used as window and nfft input
%         to spectrogram function
% 
%     whitener_rets - Structured array with fields m, j, r, c defined as
%     follows:
%         m - Variable name: mu. Double column vector of size m as
%         described in note above; contains each value of mu_k (eq 9 in
%         Helble et al (2012)) p. 2685
%         
%         j - Variable name: j_star. Double column vector of size m as
%         described in note above; contains each value of j^* (shortly
%         after eq 11 in Helble et al (2012)), p. 2685
% 
%         r - Variable name: rows. 1 x 1 Double equal to m as described in
%         note above; number of rows in X
%         
%         c - Variable name: cols. 1 x 1 Double equal to n as described in
%         note above; number of columns in X
%         
%     matrices - Structured array with fields a, b, n, nt defined as
%     follows:
%         A - Variable name: A. Double matrix of size m x n as described in
%         note above; defined by eq 7 in Helble et al (2012) p. 2685
%         
%         B - Variable name: B. Double matrix of size m x n as described in
%         note above; defined by eq 8 in Helble et al (2012) p. 2685
%         
%         N - Variable name: N. Double matrix of size m x n as described in
%         note above; defined by eq 6 in Helble et al (2012) p. 2685
%
%     X_masked - Double matrix of size m x n as describd in the note above;
%     Values are from X that lie within the bounds of a detected signal,
%     are above the threshold mu_0, and in the largest connected component
%     of the matrix adjacency graph of kept values
%
% Other Variables
%     
%     data_filtered - Double array of size samp_rate * bounds_len; result
%     of lowpass filtering data_bounded according to pass_band
%     
%     bounds_len - 1 x 1 Double, length in time of selected time series
%     specified by t_bounds, used to set spectrogram parameters to 0.05 s
%     per snapshot, 75 percent overlap
%     
%     fourier - Complex double matrix of size m' x n, as described in note
%     above; Fourier matrix output from spectrogram function
%     
%     freqs - Double array of size m' as described in note above;
%     frequencies output from spectrogram function defining the frequency
%     corresponding to each row of fourier. Values in Hz
%     
%     hz_per_bin - 1 x 1 Double, difference between adjacent values in
%     freqs, used to 'trim' freqs and fourier according to stop_band
%     
%     trim_rows - 1 x 2 Double array, row numbers of freqs/fourier
%     corresponding to stop_band     
%     
%     fourier_trimmed - Complex double matrix of size m x n, as described
%     in note above; submatrix of fourier determined by trim_rows
%     
%     c - ColorBar, value changes with each spectrogram plotting but not
%     used anywhere except for labeling the plotted ColorBar, hence its
%     variable assignment
    
function [sound, filters, original, whitener_rets, matrices, X_s, intervals, ...
    X_masked, freq_intervals_sec, noise_intervals] = GPL(fnam, wav_dir, pass_band, stop_band, ...
    t_bounds, gamma, v1, v2, eta_thresh, eta_noise, t_min, noise_thresh, ...
    max_noise_dur, min_noise_dur, filter_order, detector_type, plot_GPL_data)
    
%     Step 0: Check validity of inputs, change as necessary
    disp('Checking validity of inputs');

    [t_bounds, stop_band, pass_band, eta_thresh, eta_noise, t_min, noise_thresh, ...
        max_noise_dur, min_noise_dur] = check_inputs(fnam, wav_dir, pass_band, ...
        stop_band, t_bounds, gamma, v1, v2, eta_thresh, eta_noise, t_min, ...
        noise_thresh, max_noise_dur, min_noise_dur, filter_order, detector_type);
    
    [data, samp_rate] = audioread(strcat(wav_dir, fnam));
    
    disp('Generating Fourier matrix');
    % Step 1: Generate Fourier matrix from .wav file - Details according to
    % Helble et al (2012), p. 2690
    %     Step 1.1: Read in .wav file, assign relevant time period
    %     according to t_bounds to data
    data_bounded = data(1 + round(t_bounds(1) * samp_rate) : round(t_bounds(2) * samp_rate));

    sound.d = data_bounded;
    sound.s = samp_rate;
    %     Step 1.2: Filter data according to pass_band
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

    %     Step 1.3: Generate Fourier matrix using spectrogram function from
    %     signal processing toolbox stft using slightly different
    %     paramaters from Helble et al. (2012) p. 2690 as suggested by Dr.
    %     Spiesberger
    nfft_val = 2 ^ nextpow2(0.05 * samp_rate);
    [fourier, freqs, times] = spectrogram(data_filtered, nfft_val, 0.75 * nfft_val, ...
        nfft_val, samp_rate, 'yaxis');
         
    %     Step 1.4: Crop Fourier matrix to pass_band
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

    original.X = X;
    original.f = freqs_trimmed;
    original.t = times;
    original.n = nfft_val;
    original.h = hz_per_bin;

    disp('Determining initial test statistic');
    % Step 2: Generate matrices A, B, and N from Fourier Matrix
    %     Step 2.1: Determine conditional whitening vector mu
    [mu, j_star, rows, cols] = whitener(X);
    whitener_rets.m = mu;
    whitener_rets.j = j_star;
    whitener_rets.r = rows;
    whitener_rets.c = cols;

    %     Step 2.2: Determine matrices A, B, N
    % A is fixed-column whitening, B is fixed-row whitening
    [A, B, N] = test_stat(X, mu, gamma, v1, v2);

    matrices.A = A;
    matrices.B = B;
    matrices.N = N;
    
    disp('Running the detector procedure');
    % Step 3: Run the detector
    if detector_type == 1
        [X_s, intervals] = detector(X, gamma, v1, v2, eta_thresh, eta_noise, ...
            t_min, t_bounds, N, mu);
    else 
        intervals = detector_simple(X, eta_thresh, eta_noise, t_min, t_bounds, N);
        X_s = [];
    end
    
    disp('Running the masking procedure');
    % Step 4: Run the masking procedure
    X_masked = mask(X, intervals.i);
    
    disp('Determining frequency bounds of each signal');
    freq_intervals = box_freq(X_masked, intervals.i);
    freq_intervals_sec = (freq_intervals * hz_per_bin) + original.f(1);
    
    % Step 5: Determine noise bounds
    disp('Determining noise bounds of each signal');
    noise_intervals = noise_bounds(intervals.i, sum(N), cols, noise_thresh, ...
        t_bounds, max_noise_dur, min_noise_dur);
    
    %disp('Plotting Data');
    % Step 6: Generate Plot_Data plots
    if plot_GPL_data
        Plot_Data(original, mu, N, t_bounds, pass_band, stop_band, gamma, v1, ...
            v2, eta_thresh, eta_noise, intervals.t, X_masked, freq_intervals_sec, noise_intervals.t);
    end
end
    
    
    
    
    
    
