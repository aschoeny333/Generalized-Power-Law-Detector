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


            % Step 3: Determine the largest island of kept values in X_sig
                interval_length = length(X_sig(1, :));
                graph_edges = zeros(2, rows * interval_length);
                edge_counter = 1;
                visited_tracker = zeros(size(X_sig));
                key_dict_defs = zeros(1, numel(X_sig));
                key_counter = 1;
                for r = 1:rows
                    for c2 = 1:interval_length
                        if r > 1 && X_sig(r - 1, c2) > 0
                            if visited_tracker(r - 1, c2) == 0
                                key_dict_defs(1, key_counter) = (r - 1) * interval_length + c2;
                                visited_tracker(r - 1, c2) = key_counter;
                                key_counter = key_counter + 1;
                            end
                            if visited_tracker(r, c2) == 0
                                key_dict_defs(1, key_counter) = r * interval_length + c2;
                                visited_tracker(r, c2) = key_counter;
                                key_counter = key_counter + 1;
                            end
                            graph_edges(1, edge_counter) = visited_tracker(r, c2);
                            graph_edges(2, edge_counter) = visited_tracker(r - 1, c2);
                            edge_counter = edge_counter + 1;
                        end
                        if c2 > 1 && X_sig(r, c2 - 1) > 0
                            if visited_tracker(r, c2 - 1) == 0
                                key_dict_defs(1, key_counter) = r  * interval_length + c2 - 1;
                                visited_tracker(r, c2 - 1) = key_counter;
                                key_counter = key_counter + 1;
                            end
                            if visited_tracker(r, c2) == 0
                                key_dict_defs(1, key_counter) = r * interval_length + c2;
                                visited_tracker(r, c2) = key_counter;
                                key_counter = key_counter + 1;
                            end
                            graph_edges(1, edge_counter) = visited_tracker(r, c2);
                            graph_edges(2, edge_counter) = visited_tracker(r, c2 - 1);
                            edge_counter = edge_counter + 1;
                        end
                        if r == rows - 1 && X_sig(r + 1, c2) > 0
                            if visited_tracker(r + 1, c2) == 0
                                key_dict_defs(1, key_counter) = (r + 1) * interval_length + c2;
                                visited_tracker(r + 1, c2) = key_counter;
                                key_counter = key_counter + 1;
                            end
                            if visited_tracker(r, c2) == 0
                                key_dict_defs(1, key_counter) = r * interval_length + c2;
                                visited_tracker(r, c2) = key_counter;
                                key_counter = key_counter + 1;
                            end
                            graph_edges(1, edge_counter) = visited_tracker(r, c2);
                            graph_edges(2, edge_counter) = visited_tracker(r + 1, c2);
                            edge_counter = edge_counter + 1;
                        end
                        if c2 == interval_length - 1 && X_sig(r, c2 + 1) > 0
                            if visited_tracker(r, c2 + 1) == 0
                                key_dict_defs(1, key_counter) = r * interval_length + c2 + 1;
                                visited_tracker(r, c2 + 1) = key_counter;
                                key_counter = key_counter + 1;
                            end
                            if visited_tracker(r, c2) == 0
                                key_dict_defs(1, key_counter) = r * interval_length + c2;
                                visited_tracker(r, c2) = key_counter;
                                key_counter = key_counter + 1;
                            end
                            graph_edges(1, edge_counter) = visited_tracker(r, c2);
                            graph_edges(2, edge_counter) = visited_tracker(r, c2 + 1);
                            edge_counter = edge_counter + 1;
                        end
                    end
                end

                graph_edges(:, edge_counter:end) = [];
                key_dict_defs(:, key_counter:end) = [];

            %     Step 3.1: Determine the graph representation of kept values in X_sig
                X_sig_graph = graph(graph_edges(1, :), graph_edges(2, :));

            %     Step 3.2: Determine the largest connected component of X_sig_graph
                [bins, binsize] = conncomp(X_sig_graph);
                max_comp_ind = 1;
                max_comp_size = max(binsize);   

            %     Step 3.3: Set elements not in largest_comp to 0
                X_sig_all_bins = X_sig;
                for i = 1:length(bins)
                    if binsize(bins(i)) ~= max_comp_size
                        if mod(key_dict_defs(i), interval_length) == 0
                            X_sig(double(idivide(int16(key_dict_defs(i)), ...
                                int16(interval_length))) - 1, interval_length) = 0;
                        else
                            X_sig(double(idivide(int16(key_dict_defs(i)), ...
                                int16(interval_length))), mod(key_dict_defs(i), ...
                                interval_length)) = 0;
                        end
                    end
                end

                sig_counter = sig_counter + 1;

                for r_mask = 1 : rows
                    for c_mask = c : c + length(X_sig(1, :)) - 1
                        X_masked(r_mask, c_mask) = X_sig(r_mask, c_mask - c +  1);
                    end
                end
            end
        end
    end
end
    
    
    
    
    
    
    