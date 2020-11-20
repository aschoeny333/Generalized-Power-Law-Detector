% possible_range
% 
% Inputs: fnam, rec_j, sig_bound, programs_dir
% Output: possible range

% Determine receiver array
function [range, drift_ind, drift_date] = possible_range(fnam, rec_j, ...
    sig_bound, programs_dir)
    if (fnam(6) == 'A')
        clk_bias_path = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/Relevant Input Files/Inshore Clock Bias';
        array = 1;
    else 
        clk_bias_path = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/Relevant Input Files/Inshore Clock Bias';
        array = 2;
    end

    % Step 1: Determine clock bias
    cd(clk_bias_path)
    load('posterior_clk_bias_bnds.mat'); %#ok<*LOAD>
    bias_times = post_clk_bias.t_dn;
    biases = post_clk_bias.clk_bias_bnds;

    year = 2000 + str2double(fnam(10:11));
    month = str2double(fnam(12:13));
    day = str2double(fnam(14:15));
    hour = str2double(fnam(17:18));
    minute = str2double(fnam(19:20));
    sec = str2double(fnam(21:22));

    ref_date = [year, month, day, hour, minute, sec];

    % Find the index of the date-time corresponding with the .wav file in
    % question
    num_times = length(bias_times);
    time_difs = zeros(1, num_times);

    for i = 1:length(bias_times)
        time_difs(i) = etime(datevec(bias_times(i)), ref_date);
    end

    [~, drift_ind] = min(abs(time_difs));
    drift_date = datevec(bias_times(drift_ind));
    drift = biases(drift_ind, rec_j, array);

    if drift_ind == 0
        first_dif = etime(ref_date, datevec(bias_times(0)));
        if first_dif < 0
            if biases(1, rec_j, array) < 0
                drift = -1 * first_dif * drift_rate;
            else
                drift = first_dif * drift_rate;
            end
        end
    elseif drift_ind == num_times
        last_dif = etime(ref_date, datevec(bias_times(end)));
        if last_dif > 0
            if biases(end, rec_j, array) < 0
                drift = -1 * last_dif * drift_rate;
            else
                drift = last_dif * drift_rate;
            end
        end
    end

    % Step 2: Determine the distance between receivers
    rec_locs_path = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/Relevant Input Files';
    cd(rec_locs_path);
    load('rec_locs.mat');

    earth_rad = 6.371 * 10^6;

    if array == 1
        xy_locs = locs(2).bnds(1:5, :);
        z_locs = locs(2).z_bnds(1:5, :);
    else
        xy_locs = locs(2).bnds(6:10, :);
        z_locs = locs(2).z_bnds(6:10, :);
    end

    ref_xy_locs = xy_locs(1, :);
    ref_z_locs = z_locs(1, :);
    j_xy_locs = xy_locs(rec_j, :);
    j_z_locs = z_locs(rec_j, :);

    poss_x_difs = [ref_xy_locs(1) - j_xy_locs(1), ref_xy_locs(1) - j_xy_locs(2), ...
        ref_xy_locs(2) - j_xy_locs(1), ref_xy_locs(2) - j_xy_locs(2)];
    x_dif = max(abs(poss_x_difs)) * earth_rad;

    poss_y_difs = [ref_xy_locs(3) - j_xy_locs(3), ref_xy_locs(3) - j_xy_locs(4), ...
        ref_xy_locs(4) - j_xy_locs(3), ref_xy_locs(4) - j_xy_locs(4)];
    y_dif = max(abs(poss_y_difs)) * earth_rad;

    poss_z_difs = [ref_z_locs(1) - j_z_locs(1), ref_z_locs(1) - j_z_locs(2), ...
        ref_z_locs(2) - j_z_locs(1), ref_z_locs(2) - j_z_locs(2)];
    z_dif = max(abs(poss_z_difs));

    dist = (x_dif^2 + y_dif^2 + z_dif^2)^0.5;

    speed_bounds = [1450 1480]; % In m/s - MAKE THIS AN INPUT
    min_speed = speed_bounds(1);

    range = [sig_bound(1) - drift - dist / min_speed, sig_bound(2) - drift + dist / min_speed];

    cd(programs_dir);