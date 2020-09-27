function [new_X_sig, graph_edges, visited_tracker, key_dict, counters, ...
    X_sig_graph, bins, binsize] = main_unit(X_sig, rows)

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

counters.e = edge_counter;
counters.k = key_counter;

%  Step 3.3: Trim zeros of graph_edges and key_dict 
graph_edges(:, edge_counter:end) = [];
key_dict(:, key_counter:end) = [];

X_sig_graph = graph(graph_edges(1, :), graph_edges(2, :));

%  Step 3.4: Determine the largest connected component of X_sig_graph
[bins, binsize] = conncomp(X_sig_graph);
max_comp_size = max(binsize);   

%  Step 3.5: Set elements not in largest_comp to 0
for i = 1:length(bins)
    if binsize(bins(i)) ~= max_comp_size
        if mod(key_dict(i), interval_length) == 0
            X_sig(double(idivide(int16(key_dict(i)), ...
                int16(interval_length))) - 1, interval_length) = 0;
        else
            X_sig(double(idivide(int16(key_dict(i)), ...
                int16(interval_length))), mod(key_dict(i), ...
                interval_length)) = 0;
        end
    end
end

new_X_sig = X_sig;