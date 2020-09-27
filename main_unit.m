% Identifying and Masking the Main Spectral Content Unit of a Single
% Detected Signal
% 
% Author: Alex Schoeny
% 
% Goal: Given the submatrix defined by the time bounds of a given detected
% signal, where values below the threshold 5*mu_0 are set to 0, define a
% graph where vertices are nonzero matrix components and edges are defined
% by adjacency in the matrix. Determine the largest connected component of
% this graph and set all vertices not in it to 0
%
%     A note on size definitions: Some matrices and arrays sizes are
%     defined using m, n, and k. The values of these will change based on
%     the values of the inputs of GPL.m. Note also that m' <= m. Matrices
%     all have rows with constant frequency and columns with constant time.
% 
% Inputs
%     X_sig - Double matrix of size m x v as described in the note above,
%     where v is variable depending on inputs to the program and value of
%     sig_counter. Submatrix of X containing all rows and only columns
%     where a signal was detected according to sig_intervals
%
%     rows - 1 x 1 double equal to m as described in the note above; equal
%     to number of rows in X_sig
% 
% Outputs
%     new_X_sig - Double matrix of size m x v as described in the note
%     above, where v is variable depending on inputs to mask.m and value of
%     sig_counter. Values of the matrix are identical to X_sig, unless
%     value is less than 5*mu_0 or not in main spectral content unit
%
%     graph_edges - Double matrix of size 2 x w as described in the note
%     above, contains the endpoints of the w edges in X_sig_graph.
%     Undirected graph so no reason an endpoint is on first or second row
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
%     counters - Structured array with fields e, k defined as:
%         e - Variable name: edge_counter. 1 x 1 double looping parameter
%         tracking the number of edges (nonzero columns) in graph_edges so
%         the zeros can get trimmed off
% 
%         k - Variable name: key_counter. 1 x 1 double looping parameter
%         incrementing the number of keys in key_dict
%      
%     X_sig_graph - graph with edge_counter - 1 edges, determined from
%     adjacency of nonzero components of X_sig, largest connected component
%     determines main spectral content unit
%
%     bins - 1 x u Double matrix, output of connected components function
%     conncomp. Gives the component label of the i'th vertex in the i'th
%     position
%     
%     binsize - 1 x t Double matrix, output of connected components
%     function conncomp. Gives the size of component label i in the i'th
%     position
%
% Other Variables
%     interval_length - 1 x 1 double of size l as described in the note
%     above and in X_sig definition, used to define number of iterations in
%     part of the looping procedure
%         
%     max_comp_size - 1 x 1 Double, maximum value of binsize

function [new_X_sig, graph_edges, visited_tracker, key_dict, counters, ...
    X_sig_graph, bins, binsize] = main_unit(X_sig)

    rows = length(X_sig(:, 1));
    interval_length = length(X_sig(1, :));
    graph_edges = zeros(2, rows * interval_length);
    edge_counter = 1;
    visited_tracker = zeros(size(X_sig));
    key_dict = zeros(1, numel(X_sig));
    key_counter = 1;

    %     Step 3.2: Loop through X_sig finding edges of binary
    %     matrix and determine entries of key_dict and
    %     visited_tracker
    for r = 1:rows
        for c2 = 1:interval_length
            if X_sig(r, c2) > 0
                % Check down and to the left for all matrix
                % positions
                if r > 1 && X_sig(r - 1, c2) > 0
                    if visited_tracker(r - 1, c2) == 0
                        key_dict(1, key_counter) = (r - 1) * interval_length + c2;
                        visited_tracker(r - 1, c2) = key_counter;
                        key_counter = key_counter + 1;
                    end
                    if visited_tracker(r, c2) == 0
                        key_dict(1, key_counter) = r * interval_length + c2;
                        visited_tracker(r, c2) = key_counter;
                        key_counter = key_counter + 1;
                    end
                    graph_edges(1, edge_counter) = visited_tracker(r, c2);
                    graph_edges(2, edge_counter) = visited_tracker(r - 1, c2);
                    edge_counter = edge_counter + 1;
                end
                if c2 > 1 && X_sig(r, c2 - 1) > 0
                    if visited_tracker(r, c2 - 1) == 0
                        key_dict(1, key_counter) = r  * interval_length + c2 - 1;
                        visited_tracker(r, c2 - 1) = key_counter;
                        key_counter = key_counter + 1;
                    end
                    if visited_tracker(r, c2) == 0
                        key_dict(1, key_counter) = r * interval_length + c2;
                        visited_tracker(r, c2) = key_counter;
                        key_counter = key_counter + 1;
                    end
                    graph_edges(1, edge_counter) = visited_tracker(r, c2);
                    graph_edges(2, edge_counter) = visited_tracker(r, c2 - 1);
                    edge_counter = edge_counter + 1;
                end
                % Check above and to the right of matrix positions
                % only if directly below or to the left of the edge
                % of the matrix. This checking scheme ensures that
                % graph edges are not double counted, thus
                % requiring less pre-processing before graph
                % algorithms
                if r == rows - 1 && X_sig(r + 1, c2) > 0
                    if visited_tracker(r + 1, c2) == 0
                        key_dict(1, key_counter) = (r + 1) * interval_length + c2;
                        visited_tracker(r + 1, c2) = key_counter;
                        key_counter = key_counter + 1;
                    end
                    if visited_tracker(r, c2) == 0
                        key_dict(1, key_counter) = r * interval_length + c2;
                        visited_tracker(r, c2) = key_counter;
                        key_counter = key_counter + 1;
                    end
                    graph_edges(1, edge_counter) = visited_tracker(r, c2);
                    graph_edges(2, edge_counter) = visited_tracker(r + 1, c2);
                    edge_counter = edge_counter + 1;
                end
                if c2 == interval_length - 1 && X_sig(r, c2 + 1) > 0
                    if visited_tracker(r, c2 + 1) == 0
                        key_dict(1, key_counter) = r * interval_length + c2 + 1;
                        visited_tracker(r, c2 + 1) = key_counter;
                        key_counter = key_counter + 1;
                    end
                    if visited_tracker(r, c2) == 0
                        key_dict(1, key_counter) = r * interval_length + c2;
                        visited_tracker(r, c2) = key_counter;
                        key_counter = key_counter + 1;
                    end
                    graph_edges(1, edge_counter) = visited_tracker(r, c2);
                    graph_edges(2, edge_counter) = visited_tracker(r, c2 + 1);
                    edge_counter = edge_counter + 1;
                end
            end
        end
    end

    counters.e = edge_counter;
    counters.k = key_counter;

    %  Step 3.3: Trim zeros of graph_edges and key_dict 
    graph_edges(:, edge_counter:end) = [];
    key_dict(:, key_counter:end) = [];

    X_sig_graph = graph(graph_edges(1, :), graph_edges(2, :));

    %  Step 3.4: Determine the largest connected component of X_sig_graph
    [bins, binsize] = conncomp(X_sig_graph);
    max_comp_size = max(binsize);
    if max_comp_size < 2
        max_comp_size = -1; % May not be necessary but accounts for singleton max_comp issues
    end

    %  Step 3.5: Set elements not in largest_comp to 0 in new_X_sig
    new_X_sig = zeros(size(X_sig));
    for i = 1:length(bins)
        if binsize(bins(i)) == max_comp_size
            if mod(key_dict(i), interval_length) == 0
                r = double(idivide(int16(key_dict(i)), int16(interval_length))) - 1;
                new_X_sig(r, interval_length) = X_sig(r, interval_length);
            else
                r = double(idivide(int16(key_dict(i)), int16(interval_length)));
                c = mod(key_dict(i), interval_length);
                new_X_sig(r,c) = X_sig(r,c);
            end
        end
    end
end