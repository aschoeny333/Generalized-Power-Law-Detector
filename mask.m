% Identifying and Masking the Main Spectral Content Units of a Spectrogram
%
% Author: Alex Schoeny
% 
% Goal: Given a noise threshold above a value determined by whitener.m.
% Determine a binary matrix from values above and below the threshold and
% inside a detected signal interval. Find largest connected component of
% adjacent 1's to determine main spectral content units
%
%     A note on size definitions: Some matrices and arrays sizes are
%     defined using m, n, and k. The values of these will change based on
%     the values of the inputs of GPL.m. Note also that m' <= m. Matrices
%     all have rows with constant frequency and columns with constant time.
% 
% Inputs
%     X - Double matrix of size m x n as described in the note above;
%     absolute value of fourier_trimmed as described in GPL.m. Not
%     identical to the X defined in Helble et al. (2012) p. 2684, but the X
%     in Helble et al (2012) is only ever used in a context where its
%     absolute value is taken, hence the shorthand used here
%     
%     sig_intervals - Double matrix of size 2 x k, where each column is
%     of the form [start; end], where start and end define the beginning
%     and ending time bins of the signal identified by the detector after
%     combining and time comparison steps
%     
% Outputs
%     
%     X_masked - Double matrix of size m x n as describd in the note above;
%     Values are from X that lie within the bounds of a detected signal,
%     are above the threshold mu_0, and in the largest connected component
%     of the matrix adjacency graph of kept values
%     
% Other Variables - Used in Looping Procedure
% 
%     rows - 1 x 1 double equal to m as described in the note above; equal
%     to number of rows in X
% 
%     cols - 1 x 1 double equal to n as described in the note above; equal
%     to number of cols in X
%     
%     sig_counter = 1 x 1 double looping parameter incrementing the current
%     signal interval of interest
%     
%     X_sig - Double matrix of size m x v as described in the note above,
%     where l is variable depending on inputs to the program and value of
%     sig_counter. Submatrix of X containing all rows and only columns
%     where a signal was detected according to sig_intervals
%     
%     X_sig_row - Double matrix of size 1 x (m*v) as described in the note
%     above and in X_sig definition, used to determine binary matrix
%     threshold mu_0 from whitener.m
%     
%     mu_0 - 1 x 1 Double, binary matrix threshold parameter as determined
%     from whitener.m as described in Part B of Helble et al. (2015)
%     
%     interval_length - 1 x 1 double of size l as described in the note
%     above and in X_sig definition, used to define number of iterations in
%     part of the looping procedure
%     
%     graph_edges - Double matrix of size 2 x w as described in the note
%     above, contains the endpoints of the w edges in X_sig_graph.
%     Undirected graph so no reason an endpoint is on first or second row
%     
%     edge_counter - 1 x 1 double looping parameter tracking the number of
%     edges (nonzero columns) in graph_edges so the zeros can get trimmed
%     off
%     
%     visited tracker - Double matrix of size m x v as described in the
%     note above, tracking matrix of entries of X_sig that are zero if
%     entry is less than 5*mu_0 and increments by 1 as new entries above
%     5*mu_0 are found through looping procedure. This process also
%     generates the "keys" of key_dict
%     
%     key_dict - Double matrix 1 x u as described in the note above, stores
%     (row, column) of entry of visited_tracker with value i in its i'th
%     position
% 
%     key_counter - 1 x 1 double looping parameter incrementing the number
%     of keys in key_dict
%     
%     bins - 1 x u Double matrix, output of connected components function
%     conncomp. Gives the component label of the i'th vertex in the i'th
%     position
%     
%     binsize - 1 x t Double matrix, output of connected components
%     function conncomp. Gives the size of component label i in the i'th
%     position
%     
%     max_comp_size - 1 x 1 Double, maximum value of binsize
    
function [X_masked] = mask(X, sig_intervals)

    % Define rows, cols, sig_counter and X_masked before starting loop
    [rows, cols] = size(X);
    sig_counter = 1;
    X_masked = zeros(rows, cols);

    % Loop through columns of X and identify next signal
    for c = 1:cols
        if sig_counter <= length(sig_intervals(1, :))
            if c == sig_intervals(1, sig_counter)
            % Step 1: Determine mu_0 from X_sig
            %     Step 1.1: Determine X_sig
                X_sig = X(:, sig_intervals(1, sig_counter) : ...
                    sig_intervals(2, sig_counter));

            %     Step 1.2: Determine X_sig_row
                X_sig_row = reshape(X_sig, 1, []);

            %     Step 1.3: Determine mu_0 by whitener(X_sig_row)
                mu_0 = whitener(X_sig_row);

            % Step 2: Remove values less than 5 * mu_0 from X_sig
                for r = 1:rows
                    for c2 = 1:length(X_sig(1, :))
                        if X_sig(r, c2) < 5 * mu_0
                            X_sig(r, c2) = 0;
                        end
                    end
                end

                new_X_sig = main_unit(X_sig, rows);

                % Step 3.6: Add adjusted X_sig to X_masked
                for r_mask = 1 : rows
                    for c_mask = c : c + length(X_sig(1, :)) - 1
                        X_masked(r_mask, c_mask) = new_X_sig(r_mask, ...
                            c_mask - c +  1);
                    end
                end
                
                % Step 3.7: Increment sig_counter
                sig_counter = sig_counter + 1;
         
            end
        end
    end
end

    
    
    
    
    
    