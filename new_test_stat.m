% GPL Test Statistic Matrix Determination Function
% 
% Author: Alex Schoeny, Dr. John Spiesberger
% 
% Goal: Take the original Fourier matrix, the whitening vector determined
% from whitener.m, and the exponential parameters gamma, v1, and v2, and
% determine matrices A, B, and N according to equations 6-8 in Helble et
% al. (2012)
%
% References:
% Helble, Tyler A et al. ?A generalized power-law detection algorithm for
% humpback whale vocalizations.? The Journal of the Acoustical Society of
% America vol. 131,4 (2012): 2682-99. doi:10.1121/1.3685790
% 
%     A note on size definitions: Some matrices and arrays sizes are
%     defined using m, n, and m'. The values of these will change based on
%     the values of GPL.m, especially pass_band, stop_band, and t_bounds.
%     Matrices all have rows with constant frequency and columns with
%     constant time.
%
% Inputs: X, mu, gamma, v1, v2 from GPL.m
% 
%     X - Variable name: X. Double matrix of size m x n as described in the
%     note above; Any matrix can be used here, so m and n here may be
%     different than that of GPL.m
% 
%     mu - Variable name: mu. Double column vector of size m as described
%     in note above; contains each value of mu_k (eq 9 in Helble et al.
%     (2012)) p.2685. Must be computed from input X using whitener.m
%
%     gamma - 1 x 1 Double, paramater used to determine A and B (e.g. 1) as
%     used in Helble et al. (2012) equations 7 and 8, p.2685
% 
%     v1 - 1 x 1 Double, exponential parameter applied to A to determine N
%     (e.g. 1) as defined in Helble et al. (2012) equation 6, p.2685
% 
%     v2 - 1 x 1 Double, exponential parameter applied to B to determine N
%     (e.g. 2)as defined in Helble et al. (2012) equation 6, p.2685
% 
% Outputs:
%
%     A - Double matrix of size m x n as described in note above; defined
%     by eq 7 in Helble et al. (2012) p.2685
%         
%     B - Double matrix of size m x n as described in note above; defined
%     by eq 8 in Helble et al. (2012) p.2685
%     
%     N - Double matrix of size m x n as described in note above; defined
%     by eq 6 in Helble et al. (2012) p.2685


function [A, B, N] = new_test_stat(X, mu, gamma, v1, v2)
    X_whitened = X.^gamma - mu;
    
    % Generate matrix A
    denom_arg_A = sqrt(sum(X_whitened.^ 2, 1));
    A = abs(X_whitened) ./ denom_arg_A;
    
    % Generate matrix B
    denom_arg_B = sqrt(sum(X_whitened.^ 2, 2));
    B = abs(X_whitened) ./ denom_arg_B;
    
    % Generate test statistic
    N = (A .^ (2 * v1)) .* (B .^ (2 * v2));
