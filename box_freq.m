% Identify the frequency bounds of the masked signal
% 
% Author: Alex Schoeny
% 
% Goal: Loop through elements of X_masked to find the lowest and highest frequency bin attained by the signal. Basis for input into associator
% 
% Inputs
%     X_masked - Double matrix of size m x n as describd in the note above;
%     Values are from X that lie within the bounds of a detected signal,
%     are above the threshold mu_0, and in the largest connected component
%     of the matrix adjacency graph of kept values
%
%     sig_intervals - Double matrix of size 2 x k, where each column is
%     of the form [start; end], where start and end define the beginning
%     and ending time bins of the signal identified by the detector after
%     combining and time comparison steps
%  
% Outputs 
%     sig_intervals - Double matrix of size 2 x k, where each column is of
%     the form [start; end], where start and end define the beginning and
%     ending frequency bins of the signal contained in X_masked. Useful as
%     input to associator
% 
function [freq_intervals] = box_freq(X_masked, sig_intervals)
    rows = length(X_masked(:, 1));
    freq_intervals = zeros(size(sig_intervals));
    if ~isempty(sig_intervals)
        for i = 1:length(sig_intervals(1, :))
            min_freq = rows + 1;
            max_freq = 0;
            X_sig = X_masked(:, sig_intervals(1, i):sig_intervals(2, i));
            for r = 1:rows
                for c = 1:length(X_sig(1, :))
                    if X_sig(r,c) > 0
                        if r < min_freq
                            min_freq = r;
                        end
                        if r > max_freq
                            max_freq = r;
                        end
                    end
                end
            end
            freq_intervals(:, i) = [min_freq; max_freq];
        end
    end
end