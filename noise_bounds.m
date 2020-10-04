% Identify the noise bounds of the masked signal
% 
% Author: Alex Schoeny
% 
% Goal: Loop through elements of X_masked to find the lowest and highest frequency bin attained by the signal. Basis for input into associator
% 
% Inputs
%
%     sig_intervals - Double matrix of size 2 x k, where each column is
%     of the form [start; end], where start and end define the beginning
%     and ending time bins of the signal identified by the detector after
%     combining and time comparison steps
%  
% Outputs 
%     noise_intervals - Double matrix of size 2 x k, where each column is
%     of the form [start; end], where start and end define the beginning
%     and ending time bins of the signal and associated noise contained in
%     X_masked. Useful as input to associator
% 
function [noise_intervals] = noise_bounds(sig_intervals, test_stat, cols, ...
    noise_thresh, bins_per_sec, max_length)
    noise_intervals = zeros(size(sig_intervals));
    if ~isempty(sig_intervals)
        for i = 1:length(sig_intervals(1, :))
            dif_ind_below = 0;
            dif_ind_above = 0;
            while test_stat(sig_intervals(1, i) - dif_ind_below) > noise_thresh
                if dif_ind_below == max_length * bins_per_sec
                    break
                end
                dif_ind_below = dif_ind_below + 1;
                if sig_intervals(1, i) - dif_ind_below == 1
                    break
                end
            end
            while test_stat(sig_intervals(2,i) + dif_ind_above) > noise_thresh
                if dif_ind_above == max_length * bins_per_sec
                    break
                end
                dif_ind_above = dif_ind_above + 1;
                if sig_intervals(2, i) + dif_ind_above == cols
                    break
                end
            end
            noise_intervals(:, i) = [sig_intervals(1, i) - dif_ind_below; ...
                sig_intervals(2, i) + dif_ind_above];
        end
    end
end