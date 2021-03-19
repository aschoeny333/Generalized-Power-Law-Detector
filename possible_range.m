% Time Range Containing a Signal on a Reference Receiver
% 
% Author: Alex Schoeny
% 
% Goal: Using a detection on a reference receiver, the distance between
% receivers, the speed of sound through water, clock bias measurements to
% synchronize receivers, determine the time bounds on another receiver
% which must contain a detection from a reference receiver.
% 
% Inputs: fnam, rec_j, sig_bound, programs_dir
%
%     fnam - Char array, .wav file name (e.g. 'test.wav')
%
%     rec_j - 1 x 1 double, integer value for which receiver to correlate
%     with
%
%     sig_bound - 1 x 2 double, time bounds of the relevant detection on
%     the reference receiver (in sec)
%
%     programs_dir - Char array, gives directory containing associator.m
%     and other relevant programs
%
% Output: 
%
%     range - 1 x 2 Double, time bounds on non-reference receiver
%     which the detection in question from reference receiver must occur.
%
%     drift_ind - 1 x 1 Double, index of clock biases array where the
%     corresponding date-time is closest to the date-time of the relevant
%     .wav file
%     
%     drift_date - 1 x 6 Double, vector containing the dete-time components
%     of the closest clock bias measurement
%
%     drift - 1 x 1 Double, value in seconds by which the time series on
%     receiver j needs to be translated to align properly with the
%     reference receiver
%     
%     dist - 1 x 1 Double. Given the possible bounds in each dimension
%     (longitude, latitude, altitude), find the maximum possible distance
%     between the two receivers in meters

function [range, drift_ind, drift_date, drift, dist] = possible_range(fnam, ...
    rec_j, sig_bound, wav_dir)
    
    wav_dir_path_length = length(wav_dir);
    % Step 1: Determine the receiver array from fnam
    if (fnam(wav_dir_path_length + 6) == 'A')
        clk_bias_path = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/Relevant Input Files/Inshore Clock Bias/';
        array = 1;
    else 
        clk_bias_path = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/Relevant Input Files/Midshore Clock Bias/';
        array = 2;
    end

    % Step 2: Determine clock bias
    % Step 2.1: Load file containing biases time series
    load(strcat(clk_bias_path, 'posterior_clk_bias_bnds.mat')); %#ok<*LOAD>
    bias_times = post_clk_bias.t_dn;
    biases = post_clk_bias.clk_bias_bnds;

    % Step 2.2: Determine the timestamp (ref_date) of the file from fnam
    year = 2000 + str2double(fnam(wav_dir_path_length + 10:wav_dir_path_length + 11));
    month = str2double(fnam(wav_dir_path_length + 12:wav_dir_path_length + 13));
    day = str2double(fnam(wav_dir_path_length + 14:wav_dir_path_length + 15));
    hour = str2double(fnam(wav_dir_path_length + 17:wav_dir_path_length + 18));
    minute = str2double(fnam(wav_dir_path_length + 19:wav_dir_path_length + 20));
    sec = str2double(fnam(wav_dir_path_length + 21:wav_dir_path_length + 22));

    ref_date = [year, month, day, hour, minute, sec];

    % Step 2.3: Find the index of the closest date-time with a recorded
    % clock bias measurement corresponding with the .wav file in question
    num_times = length(bias_times);
    time_difs = zeros(1, num_times);

    for i = 1:length(bias_times)
        time_difs(i) = etime(datevec(bias_times(i)), ref_date);
    end

    [~, drift_ind] = min(abs(time_difs));
    drift_date = datevec(bias_times(drift_ind));

    % Step 2.4: If nearest date-time corresponds to first or last recorded
    % clock bias, determine drift by linear extrapolation using drift_rate.
    % Otherwise, determine drift by linear interpolation of the two neares
    % clock bias bounds measurements
    if drift_ind == 0
        % Extrapolation procedure - reference date before any measurements
        first_dif = etime(ref_date, datevec(bias_times(0)));
        if first_dif < 0
            if biases(1, rec_j, array) < 0
                drift = -1 * first_dif * drift_rate;
            else
                drift = first_dif * drift_rate;
            end
        end
    elseif drift_ind == num_times
        % Extrapolation procedure - reference date after all measurements
        last_dif = etime(ref_date, datevec(bias_times(end)));
        if last_dif > 0
            if biases(end, rec_j, array) < 0
                drift = -1 * last_dif * drift_rate;
            else
                drift = last_dif * drift_rate;
            end
        end
    else
        % Interpolation procedure
        dates_dif = etime(drift_date, ref_date);
        if dates_dif > 0
            gap_secs = etime(drift_date, datevec(bias_times(drift_ind - 1)));
            drift = interp1([0, gap_secs], [biases(drift_ind - 1, rec_j, array), ...
                biases(drift_ind, rec_j, array)], gap_secs - dates_dif);
        elseif dates_dif < 0
            gap_secs = etime(datevec(bias_times(drift_ind + 1)), drift_date);
            drift = interp1([0, gap_secs], [biases(drift_ind, rec_j, array), ...
                biases(drift_ind + 1, rec_j, array)], -1 * dates_dif);
        else
            drift = biases(drift_ind, rec_j, array);
        end  
    end

    % Step 3: Determine the distance between receivers
    % Step 3.1: Load coordinate locations of both receivers
    load('/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/Relevant Input Files/rec_locs.mat');

    if array == 1
        xy_locs = locs(2).bnds(1:5, :);
    else
        xy_locs = locs(2).bnds(6:10, :);
    end

    ref_xy_locs = xy_locs(1, :);
    j_xy_locs = xy_locs(rec_j, :);
    
    % Step 3.2: Determine the maximum possible distance between the
    % receivers in each dimension
    poss_x_difs = [ref_xy_locs(1) - j_xy_locs(1), ref_xy_locs(1) - j_xy_locs(2), ...
        ref_xy_locs(2) - j_xy_locs(1), ref_xy_locs(2) - j_xy_locs(2)];
    [~, max_lon_dif_ind] = max(abs(poss_x_difs));
    if max_lon_dif_ind == 1
        lon1 = ref_xy_locs(1);
        lon2 = j_xy_locs(1);
    elseif max_lon_dif_ind == 2
        lon1 = ref_xy_locs(1);
        lon2 = j_xy_locs(2);
    elseif max_lon_dif_ind == 3
        lon1 = ref_xy_locs(2);
        lon2 = j_xy_locs(1);
    else
        lon1 = ref_xy_locs(2);
        lon2 = j_xy_locs(2);
    end

    poss_y_difs = [ref_xy_locs(3) - j_xy_locs(3), ref_xy_locs(3) - j_xy_locs(4), ...
        ref_xy_locs(4) - j_xy_locs(3), ref_xy_locs(4) - j_xy_locs(4)];
    [~, max_lat_dif_ind] = max(abs(poss_y_difs));
    if max_lat_dif_ind == 1
        lat1 = ref_xy_locs(3);
        lat2 = j_xy_locs(3);
    elseif max_lat_dif_ind == 2
        lat1 = ref_xy_locs(3);
        lat2 = j_xy_locs(4);
    elseif max_lat_dif_ind == 3
        lat1 = ref_xy_locs(4);
        lat2 = j_xy_locs(3);
    else 
        lat1 = ref_xy_locs(4);
        lat2 = j_xy_locs(4);
    end


    % Step 3.3: Determine the maximum Euclidean distance between the two
    % receivers
    dist = great_circle_distance(lat1, lat2, lon1, lon2);

    % Step 4: Determine the possible range on receiver j that must contain
    % the reference signal
    speed_bounds = [1450 1480]; % In m/s
    min_speed = speed_bounds(1);

    % Note: A negative drift value indicates that a receiver is slow
    % compared to the reference receiver. Thus, to account for clock drift
    % one must add that negative value to the association range to check
    % earlier values in the time series. The same thinking applies for
    % positive drift values, but in reverse.
    range = [sig_bound(1) + drift - dist / min_speed, sig_bound(2) + drift + dist / min_speed];


