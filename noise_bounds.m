% Identify the noise bounds of the masked signal
% 
% Author: Alex Schoeny
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
    noise_thresh, t_bounds, max_noise_dur)
    bins_per_sec = cols / (t_bounds(2) - t_bounds(1));
    bin_intervals = zeros(size(sig_intervals));
    if ~isempty(sig_intervals)
        for i = 1:length(sig_intervals(1, :))
            dif_ind_below = 0;
            dif_ind_above = 0;
            while test_stat(sig_intervals(1, i) - dif_ind_below) < noise_thresh
                if dif_ind_below > max_noise_dur * bins_per_sec
                    disp("Max length reached");
                    break
                end
                dif_ind_below = dif_ind_below + 1;
                if sig_intervals(1, i) - dif_ind_below == 1
                    disp("End of matrix reached");
                    break
                end
            end
            while test_stat(sig_intervals(2,i) + dif_ind_above) < noise_thresh
                if dif_ind_above > max_noise_dur * bins_per_sec
                    disp("Max length reached");
                    break
                end
                dif_ind_above = dif_ind_above + 1;
                if sig_intervals(2, i) + dif_ind_above == cols
                    disp("End of matrix reached");
                    break
                end
            end
            bin_intervals(:, i) = [sig_intervals(1, i) - dif_ind_below; ...
                sig_intervals(2, i) + dif_ind_above];
        end
    end
    
    noise_intervals.b = bin_intervals;
    noise_intervals.t = bin_intervals / bins_per_sec + t_bounds(1);
end