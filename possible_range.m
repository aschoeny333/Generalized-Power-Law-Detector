% possible_range
% 
% Inputs: Distance between receivers, signal bounds on reference receiver, clock bias, min speed of sound
% 
% Output: possible range

% Determine receiver array
if (fnam(6) == 'A')
    clk_bias_path = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/Relevant Input Files/Inshore Clock Bias';
    array = 1;
else 
    clk_bias_path = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/Relevant Input Files/Inshore Clock Bias';
    array = 2;
end

% Step 1: Determine clock bias
cd(clk_bias_path)
load('posterior_clk_bias_bnds.mat');
bias_times = post_clk_bias.t_dn;
biases = post_clk_bias.clk_bias_bnds;

year = 2000 + str2double(fnam(10:11));
month = str2double(fnam(11:12));
day = str2double(fnam(13:14));
hour = str2double(fnam(16:17));
min = str2double(fnam(18:19));
sec = str2double(fnam(20:21));

% Find the index of the date-time corresponding with the .wav file in
% question
time_index = 0;
for i = 1:length(bias_times)
    cur_date = datevec(bias_times(i));
    if day == cur_date(3)
        if hour == cur_date(4)
            if month == cur_date(2)
                if year == cur_date(1)
                    time_index = i;
                    break;
                end
            end
        end
    end
end

% If a clock bias measurement exists, assign it to drift. If not, determine
% nearest bias measurement and assign to drift its value plus additional
% possible drift according to drift_rate
if i ~= 0
    drift = biases(i, rec, array);
else
    drift_rate = 2.5*10^-6;
    first_date = datevec(bias_times(1));
    ref_date = [year, month, day, hour, min, sec];
    first_dif = etime(first_date, ref_date);
    if first_dif > 0
        if biases(1, rec, array) < 0
            drift = -1 * first_dif * drift_rate;
        else
            drift = first_dif * drift_rate;
        end
    end
    
    last_date = datevec(bias_times(end));
    last_dif = etime(last_date, ref_date);
    if last_dif < 0
        if biases(end, rec, array) < 0
            drift = -1 * first_dif * drift_rate;
        else
            drift = first_dif * drift_rate;
        end   
    else
        disp('Error: possible range procedure failed, must de-bug');
    end
end

% Step 2: Determine the distance between receivers
rec_locs_path = '/Users/Alex_Schoeny/Desktop/Research/GPL/Programs and Test Files - Dev/Relevant Input Files';
cd(rec_locs_path);
load('rec_locs.mat');

dist; % TO DO: Write script to determine dist between two receivers

speed_bounds = [1450 1480]; % In m/s
min_speed = speed_bounds(1);



range = [sig_bound(1) - drift - dist / min_speed, sig_bound(2) - drift + dist / min_speed];