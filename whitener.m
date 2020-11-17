% Conditional whitening vector determination function 
%
% Author: Alex Schoeny
%
% Goal: Determine conditional whitening vector mu as defined by Helble et
% al, as well as j_star (a vector containing the value of j^* as defined in
% Helble et al for each row), and integers rows and cols storing the size
% dimensions of the input matrix X. mu, rows, and cols are used later in
% GPL.m, while j_star may be relevant for testing purposes. 
%
% References:
% Helble, Tyler A et al. ?A generalized power-law detection algorithm for
% humpback whale vocalizations.? The Journal of the Acoustical Society of
% America vol. 131,4 (2012): 2682-99. doi:10.1121/1.3685790
% 
% Inputs: X from GPL.m
%
%     A note on size definitions: Some matrices and arrays sizes are
%     defined using m, n, and m'. The values of these will change based on
%     the values of the inputs of GPL.m, especially pass_band, stop_band,
%     and t_bounds. Matrices all have rows with constant frequency and
%     columns with constant time.
% 
%     X - Variable name: X. Double matrix of size m x n as described in the
%     note above; Any matrix can be used here, so m and n here may be
%     different than that of GPL.m
% 
% Outputs:
% 
%     mu - Double column vector of size m as described in note above;
%     contains each value of mu_k (eq 9 in Helble et al (2012)) p. 2685
%
%     j_star - Double column vector of size m' as described in note above;
%     contains each value of j^* (shortly after eq 11 in Helble et al
%     (2012)) p.2685
% 
%     rows - 1 x 1 Double equal to m' as described in note above; number of
%     rows in X
%     
%     cols - 1 x 1 Double equal to n as described in note above; number of
%     columns in X

function [mu, j_star, rows, cols] = whitener(X)
    [rows, cols] = size(X);
    % The expression below is for J/2 - 1; unclear in Helble et al if J/2
    % is rounded up or down; rounding down here
    dif_ind = fix(cols / 2) - 1; 
    S = sort(X, 2);
    % Initializing j_star and mu to be 'filled in' during matix iteration
    j_star = zeros(rows, 1);
    mu = zeros(rows, 1);
    % Iterating through matrix rows
    for r = 1:rows
        min_j = 1;
        min_dif = S(r, 1 + dif_ind) - S(r, 1);
        % Iterating through sorted matrix columns
        for c = 1:cols - dif_ind
            if S(r, c + dif_ind) - S(r, c) < min_dif
                min_dif = S(r, c + dif_ind) - S(r, c);
                min_j = c;
            end
        end
        j_star(r, 1) = min_j;
        mu(r,1) = 2 * sum(S(r, min_j : min_j + dif_ind)) / cols;
    end
    
    
