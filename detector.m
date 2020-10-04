% Generalized Power Law Detector Implementation 
% Author: Alex Schoeny 
% Goal: Using manipulated spectrogram from GPL.m, identify time intervals
% where a signal is detected 
% 
% References:
% Helble, Tyler A et al. ?A generalized power-law detection algorithm for
% humpback whale vocalizations.? The Journal of the Acoustical Society of
% America vol. 131,4 (2012): 2682-99. doi:10.1121/1.3685790
% 
% Inputs:
%
%     A note on size definitions: Some matrices and arrays sizes are
%     defined using m, n, and n'. The values of these will change based on
%     the values of the inputs of GPL.m, especially pass_band, stop_band,
%     and t_bounds. Matrices all have rows with constant frequency and
%     columns with constant time. Note that in the below descriptions, k <=
%     k' <= k'', since each step may either combine or remove one of the
%     time interval bounds. Similarly, n' <= n
%
%     X - Double matrix of size m x n as described in note above; absolute
%     value of fourier_trimmed. Not identical to the X defined in Helble et
%     al (2012) p.2684 , but the X in Helble et al (2012) is only ever used
%     in a context where its absolute value is taken, hence the shorthand
%     used here
%
%     mu - Double column vector of size m as described in note above;
%     contains each value of mu_k (eq 9 in Helble et al (2012)) p.2685
%
%     N - Double matrix of size m x n as described in note above; defined
%     by eq 6 in Helble et al (2012) p.2685
% 
%     eta_thresh - 1 x 1 Double, minimum test statistic threshold a signal
%     must reach in at least one time bin in order to be identified by
%     detector.m. Derived largely in sections IIIB-D in Helble et al. Exact
%     value given in Section IV (2.62 * 10^-4) p.2690
% 
%     eta_noise - 1 x 1 Double, minimum test statistic threshold a signal
%     must reach in all time bins in order to have its time interval
%     determined by detector.m; that is, test statistics below this value
%     imply a signal does not exist in that time bin. Derived largely in
%     sections IIIB-D in Helble et al. Exact value given in Section IV
%     (2.07 * 10^-5) p.2690
%
%     t_min - 1 x 1 Double, minimum time length of a signal in seconds
%     (e.g. 0.35)
% 
% Outputs
% 
%     X_s - Double matrix of size m x n', consisting of columns of X
%     satisfying the criteria laid out in Helble et al. (2012) page 2691
%     
%     intervals - Structured array with fields i, u, w, t defined as
%     follows:
%         u - Variable name: signal_intervals_uncombined. Double matrix of
%         size 2 x k, where each column is of the form [start; end], where
%         start and end define the beginning and ending time bins of the
%         signal identified by the detector before combining and time
%         comparison steps
%     
%         w - Variable name: signal_intervals_with_short. Double Matrix of
%         size 2 x k' matrix, where each column is of the form [start;
%         end], where start and end define the beginning and ending time
%         bins of the signal identified by the detector after combining and
%         before time comparison steps
%
%         i - Variable name: signal_intervals. Double matrix of size 2 x
%         k'', where each column is of the form [start; end], where start
%         and end define the beginning and ending time bins of the signal
%         identified by the detector after combining and time comparison
%         steps
%    
%         t - Variable name: time_intervals. Double matrix of size 2 x k'',
%         where each column is of the form [start; end], where start and
%         end define the beginning and ending time value in seconds of the
%         signal identified by the detector after combining and time
%         comparison steps
% 
% 
% Other Variables
% 
%     col_counter - 1 x 1 Double, looping variable used in multiple places
%     to count the number of columns identifed to be removed after the
%     loop's termination. This technique (filling a zeros array partially
%     then truncating the zeros after a looping procedure) is used to
%     prevent resizing, which is a strain on the program's runtime.
% 
%     mu_s - Double column vector of size m as described in note above,
%     output of whitener.m when X_s used as input, necessary for finding
%     columns of X with test statistic below eta_noise while accounting for
%     influence of signals on standard mu calculation
% 
%     N_s - Double array ending as size m*n' as described in note above,
%     test statistic matrix calculated using X_s and mu_s for purpose
%     described in mu_s description
% 
%     test_stat_s - Double row vector of size n' as described in note
%     above, test statistic array calculated from N_s for purpose described
%     in mu_s description
% 
%     remove_cols - Double row vector initially of size n as described in
%     note above but later trimmed to only contain column indeces of X_s
%     that had test statistic test_stat_s less than or equal to eta_thresh
% 
%     cols_s - 1 x 1 Double equal to n' as described in note above, length
%     of X_s used to generated a random integer to randomly select a column
%     in X_s to replace columns in X where a signal has been identified and
%     recorded
% 
%     signal_counter - 1 x 1 Double, looping variable counting number of
%     identified signals to truncate zeros after loop terminates
% 
%     cur_N - Double matrix of size m x n as described in note above,
%     looping variable that updates with every pass of the loop with the
%     new test statistic matrix calculated from cur_matrix and cur_mu
% 
%     cur_test_stat - Double row vector of size n as described in note
%     above, looping variable that updates with every pass of the loop with
%     the new test statistic calculated from cur_N
% 
%     cur_matrix - Double matrix of size m x n as described in note above
%     that updates with every pass of the loop with the new matrix after
%     swapping identified signal columns with randomly selected columns
%     from X_s
% 
%     cur_mu - Double column vector of size m as described in note above
%     that updates with every pass of the loop with the new whitener array
%     calculated from cur_matrix
% 
%     dif_ind_below - 1 x 1 Double, looping variable that starts at max_ind
%     and increments down by 1 until a test statistic value below eta_noise
%     is reached
% 
%     dif_ind_above - 1 x 1 Double, looping variable that starts at max_ind
%     and increments up by 1 until a test statistic value below eta_noise
%     is reached
% 
%     max_ts - 1 x 1 Double, maximum test statistic identified in
%     cur_test_stat. A signal will be identified if max_ts is greater than
%     eta_thresh, and the loop terminates if not
% 
%     max_ind - 1 x 1 Double, index of max_ts as desribed above
% 
%     interval_ind_below - 1 x 1 Double, index determined by subtracting
%     dif_ind_below from max_ind, determines lower bound of the signal and
%     is stored in signal_intervals
% 
%     interval_ind_above - 1 x 1 Double, index determined by adding
%     dif_ind_above to max_ind, determines upper bound of the signal and is
%     stored in signal_intervals
% 
%     sec_per_bin - 1 x 1 Double, calculated value of number of seconds in
%     each time bin, used to calculate time_intervals from signal_intervals
% 
%     min_bins - 1 x 1 Double, minimum number of time bins that must be
%     contained in a signal interval for the signal to have a duration of
%     at least t_min
% 
%     remove_intervals - Double row vector intially of length k' as
%     described in note above but later trimmed to only contain indeces of
%     signal_intervals that have duration at least t_min

function [X_s, intervals] = detector(X, gamma, v1, v2, eta_thresh, ...
    eta_noise, t_min, t_bounds)

    % Step 1: Preprocessing to determine X_s as defined on page 2691

    % Step 1.1 - Determine mu and N
    mu = whitener(X);
    [~, ~, N] = test_stat(X, mu, gamma, v1, v2);

    % Step 1.1 - Identify columns of X that have a test statistic less than
    % eta_noise
    N_sum = sum(N);
    [rows, cols] = size(N);
    X_s = zeros(rows, cols);
    col_counter = 1;
    for i = 1:cols
        if N_sum(i) < eta_noise
            X_s(:, col_counter) = X(:, i);
            col_counter = col_counter + 1;
        end
    end  

    % Step 1.2 - Determine test statistic array for X_s after removing 0
    % columns
    X_s(:, col_counter : end) = [];
    mu_s = whitener(X_s);
    [~, ~, N_s] = test_stat(X_s, mu_s, gamma, v1, v2);
    test_stat_s = sum(N_s);

    % Step 1.3 - Remove columns of X_s with test statistic greater than
    % eta_thresh
    remove_cols = zeros(cols, 1);
    col_counter = 1;
    for i = 1:length(N_s)
        if test_stat_s(i) > eta_thresh
            remove_cols(col_counter) = i;
        end
    end
    remove_cols(col_counter : end) = [];
    X_s(:, remove_cols) = [];

    % Step 2: Iterative determination of event intervals

    cols_s = length(X_s(:, 1));
    if cols_s == 0
        % What to do here? is this case possible?
        disp("X_s is empty. Unknown how to handle this case.");
    else
        % Declare return variable 
        signal_intervals = zeros(2, cols);
        signal_counter = 1;

        % Declare looping variables
        cur_test_stat = N_sum;
        cur_N = N;
        cur_matrix = X;
        cur_mu = mu;
        dif_ind_below = 0;
        dif_ind_above = 0;

        % Find current maximum test statistic value and index
        [max_ts, max_ind] = max(cur_test_stat);

        % Looping procedure: Identify all signal time segments satisfying
        % criteia laid out in Helble et al. (2012) p.2691
        while max_ts > eta_thresh
            % Display alert that the loop has restarted. Useful for
            % understanding timing of program execution
            disp("another max found");

            % Identify time bounds of signal 
            while cur_test_stat(max_ind - dif_ind_below) > eta_noise
                dif_ind_below = dif_ind_below + 1;
                if max_ind - dif_ind_below == 1
                    break
                end
            end
            while cur_test_stat(max_ind + dif_ind_above) > eta_noise
                dif_ind_above = dif_ind_above + 1;
                if max_ind + dif_ind_above == cols
                    break
                end
            end

            % Assign signal bound indeces based on dif_ind variables
            interval_ind_below = max_ind - dif_ind_below;
            interval_ind_above = max_ind + dif_ind_above;

            % Record time interval
            signal_intervals(:, signal_counter) = [interval_ind_below ; interval_ind_above];
            signal_counter = signal_counter + 1;

            % Add random X_s columns in place of signal columns 
            for i = interval_ind_below : interval_ind_above
                cur_matrix(:, i) = X_s(:, randi(cols_s));
            end

            % Reset loop variable values
            dif_ind_below = 0;
            dif_ind_above = 0;
            cur_mu = whitener(cur_matrix);
            [~, ~, cur_N] = test_stat(cur_matrix, cur_mu, gamma, v1, v2);
            cur_test_stat = sum(cur_N);
            [max_ts, max_ind] = max(cur_test_stat);
        end
    end

    % Step 3: Trim 0's from signal_intervals, combine adjacent events, check
    % intervals against t_min

    % Step 3.1: Trim 0's, sort, and save a copy for debugging
    signal_intervals(:, signal_counter : end) = [];
    signal_intervals = sort(signal_intervals, 2);
    signal_intervals_uncombined = signal_intervals;
    
    % Step 3.2: Combine adjacent events (within 1 time bin)
    for i = 1 : signal_counter - 2
        % while loop used to combine multiple adjacent events
        if i >= length(signal_intervals)
            break
        end
        while signal_intervals(1, i + 1) - signal_intervals(2, i) < 3
            signal_intervals(2, i) = signal_intervals(2, i + 1);
            signal_intervals(:, i + 1) = [];
            if i >= length(signal_intervals(1,:))
                break
            end
        end
    end
    
    % Step 3.3: Remove intervals shorter than t_min after saving copy for
    % debugging
    signal_intervals_with_short = signal_intervals;
    
    sec_per_bin = (t_bounds(2) - t_bounds(1)) / cols;
    min_bins = t_min / sec_per_bin;
    remove_intervals = zeros(length(signal_intervals), 1);
    col_counter = 1;
    for i = 1 : length(signal_intervals(1, :))
        if signal_intervals(2, i) - signal_intervals(1, i) < min_bins
            remove_intervals(col_counter) = i;
            col_counter = col_counter + 1;
        end
    end
    remove_intervals(col_counter : end) = [];
    signal_intervals(:, remove_intervals) = [];
     
    intervals.u = signal_intervals_uncombined;
    intervals.w = signal_intervals_with_short;
    intervals.i = signal_intervals;

    % Step 4: Translate intervals from bin number to seconds
    time_intervals = signal_intervals * sec_per_bin + t_bounds(1);
    intervals.t = time_intervals;







