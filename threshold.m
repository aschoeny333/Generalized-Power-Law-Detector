% Post-GPL minimum power thresholding function
% Author: Alex Schoeny
% Goal: Introduce minimum power threshold to matrix N resulting from GPL.m
% to make the background noise appear to be a more uniform color on final
% spectrogram
% Inputs: N, thresh from GPL.m
%     N - Double matrix of size m' x n, where m' is the number of frequency
%     bins (namely, trim_rows(2)-trim_rows(1)) and n is the number of
%     snapshots; defined by eq 6 in Helble et al (2012)
% 
%     thresh - 1 x 1 Double, minimum power threshold of final spectrogram.
%     Not to be confused with eta_thresh of Helble et al (2012) (e.g.
%     10^-7)
%

function N_thresh = threshold(N, thresh)
    [rows, cols] = size(N);
    N_thresh = zeros(rows, cols);
    for r = 1:rows
        for c = 1:cols
            if N(r, c) < thresh
                N_thresh(r, c) = thresh;
            else
                N_thresh(r, c) = N(r, c);
            end
        end
    end
    