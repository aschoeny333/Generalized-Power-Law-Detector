% Identify the noise bounds of the masked signal
% 
% Author: Alex Schoeny, Dr. John Spiesberger
% 
% Goal: Loop through signals identified by detector and increment away from
% either side of the signal bounds to determine noise bounds. Incrementing
% stops when noise threshold in test statistic or maximum noise length is
% reached
% 
% Inputs
%
%     sig_intervals - Double matrix of size 2 x k, where each column is
%     of the form [start; end], where start and end define the beginning
%     and ending time bins of the signal identified by the detector after
%     combining and time comparison steps
%
%     test_stat - Double row vector of size n as described in note above;
%     likely test statistic for original spectrogram defined by eq 2 in
%     Helble et al (2012) p.2684, where v = 1
%
%     cols - Variable name: cols. 1 x 1 Double equal to n as described in
%     note above; number of columns in X
%
%     noise_thresh - 1 x 1 double, value of test statistic above which
%     noise interval iteration procedure stops. Not referenced in either
%     Helble paper, but consider using eta_noise
%
%     t_bounds - 1 x 2 Double array, time interval of interest of form
%     [start time, end time] in sec (e.g. [360, 450]). If start time > end
%     time, values are swapped. If any other improper entry (e.g. start
%     time = end time, start time < 0, end time > time series length),
%     t_bounds is set to [0, time series length]
% 
%     max_noise_dur - 1 x 1 double, value in seconds that the length of the
%     noise bounds cannot exceed
%  
% Outputs 
%     noise_intervals - Structured array with fields b, t defined as
%     follows:
%         b - Variable name: bin_intervals. Double matrix of size 2 x k,
%         where each column is of the form [start; end], where start and
%         end define the beginning and ending time bins (column values) of
%         the signal and associated noise contained in X_masked. Useful as
%         input to associator 
%
%         t - No variable name, defined directly from bin_intervals. Double
%         matrix of size 2 x k, where each column is of the form [start;
%         end], where start and end define the beginning and ending time in
%         seconds of the signal and associated noise contained in X_masked.
%         Useful as input to associator and Plot_Data.m
    
function [noise_intervals] = noise_bounds(sig_intervals, test_stat, cols, ...
    noise_thresh, t_bounds, max_noise_dur, min_noise_dur)

    bins_per_sec = cols / (t_bounds(2) - t_bounds(1));
    bin_intervals = zeros(size(sig_intervals));
    num_sigs = length(sig_intervals(1, :));
    
    if ~isempty(sig_intervals)
        for i = 1:num_sigs
            
            % Define iteration variables
            dif_ind_below = 0;
            dif_ind_above = 0;
            
            % Determine if there is a previous signal and if so assign
            % prev_sig to its bin value. If not, assign it to the beginning
            % of the matrix
            if i - 1 > 0
                prev_sig = sig_intervals(2, i - 1);
            else 
                prev_sig = 1;
            end
            
            % Determine if there is a next signal and if so assign next_sig
            % to its bin value. If not, assign it to the end of the matrix
            if i + 1 <= num_sigs
                next_sig = sig_intervals(1, i + 1);
            else 
                next_sig = cols;
            end
            
            % Iterate backwards from the start of the signal
            while test_stat(sig_intervals(1, i) - dif_ind_below) < noise_thresh
                % Check if noise bounds are too long yet
                if dif_ind_below > max_noise_dur * bins_per_sec
                    disp('Max length reached');
                    break
                end
                
                % Check if noise bounds encroaching on previous signal or
                % beginning of matrix reached
                if sig_intervals(1, i) - dif_ind_below == prev_sig
                    disp('Previous signal reached');
                    break
                end
                
                if sig_intervals(1, i) - dif_ind_below == 1
                    disp('Beginning of the audio file reached');
                    break
                end
                
                % Increment dif_ind_below to continue iteration
                dif_ind_below = dif_ind_below + 1;
            end
            
            % Iterate forwards from the end of the signal
            while test_stat(sig_intervals(2,i) + dif_ind_above) < noise_thresh
                % Check if noise bounds are too long yet
                if dif_ind_above > max_noise_dur * bins_per_sec
                    disp('Max length reached');
                    break
                end
                
                % Check if noise bounds encroaching on next signal or end
                % of matrix reached
                if sig_intervals(2, i) + dif_ind_above == next_sig
                    disp('Next signal reached');
                    break
                end
                
                if sig_intervals(2, i) + dif_ind_above == length(test_stat)
                    disp('End of audio file reached');
                    break
                end
         
                % Increment dif_ind_above to continue iteration
                dif_ind_above = dif_ind_above + 1;
            end
            % Case: If duration of combined noise regions is less than
            % min_noise_dur, determine remaining duration, iterate through
            % b + 1 possible selections (where b is bins equiv to necessary
            % duration) of noise region, determine min f(b-region_i), set
            % as new bin_intervals entry
            
            % DISCUSS: Optimum function to evaluate optimum x_test (sum,
            % max, etc.)
            
            if (dif_ind_below + dif_ind_above) < min_noise_dur * bins_per_sec
                % CHECK some of these might have decimals, screwing things
                % up?
                dur_remaining = round(min_noise_dur * bins_per_sec - ...
                    (dif_ind_below + dif_ind_above));
                optimum_bound = 0;
                optimum_bound_val = 10^5;
                for x = 0:dur_remaining
                    x_test = 0;
                    for y = 0:(dur_remaining - x)
                        x_test = x_test + test_stat(sig_intervals(1, i) - ...
                            dif_ind_below - dur_remaining + x + y);
                    end
                    if x > 0
                        for y = 1:x
                            x_test = x_test + test_stat(sig_intervals(2, i) + ...
                            dif_ind_above + y);
                        end
                    end
                    if x_test < optimum_bound_val
                        optimum_bound = x;
                    end
                end
                
                % Append optimally extended noise bounds to return variable 
                bin_intervals(:, i) = [sig_intervals(1, i) - dif_ind_below - ...
                    dur_remaining + optimum_bound; sig_intervals(2, i) + ...
                    dif_ind_above + optimum_bound];
            else
                % Append noise bounds to return variable
                bin_intervals(:, i) = [sig_intervals(1, i) - dif_ind_below; ...
                    sig_intervals(2, i) + dif_ind_above];
            end

        end
    end
    
    % Assign fields of return variable
    noise_intervals.b = bin_intervals;
    noise_intervals.t = bin_intervals / bins_per_sec + t_bounds(1);
end