test_empty = zeros(4);
tested_empty = main_unit(test_empty);
assert(isequal(tested_empty, zeros(4)));

test_no_adj = [[0, 1, 0, 1]; [1, 0, 1, 0]; [0, 1, 0, 1]; [1, 0, 1, 0]];
[tested_no_adj, adj_edges, ~, ~, ~, empty_graph] = main_unit(test_no_adj);
assert(isequal(tested_no_adj, zeros(4)));

test_interior = [[1, 0, 0, 0]; [0, 1, 1, 0]; [0, 1, 1, 0]; [0, 0, 0, 1]];
test_interior_check = [[0, 0, 0, 0]; [0, 1, 1, 0]; [0, 1, 1, 0]; [0, 0, 0, 0]];
tested_interior = main_unit(test_interior);
assert(isequal(tested_interior, test_interior_check));

test_top_edge = [[1, 0, 0, 0]; [0, 1, 1, 0]; [0, 1, 1, 0]; [0, 1, 1, 0]];
test_top_edge_check = [[0, 0, 0, 0]; [0, 1, 1, 0]; [0, 1, 1, 0]; [0, 1, 1, 0]];
tested_top_edge = main_unit(test_top_edge);
assert(isequal(tested_top_edge, test_top_edge_check));

test_bottom_edge = [[0, 1, 1, 0]; [0, 1, 1, 0]; [0, 1, 1, 0]; [0, 0, 0, 1]];
test_bottom_edge_check = [[0, 1, 1, 0]; [0, 1, 1, 0]; [0, 1, 1, 0]; [0, 0, 0, 0]];
tested_bottom_edge = main_unit(test_bottom_edge);
assert(isequal(tested_bottom_edge, test_bottom_edge_check));

test_left_edge = [[0, 0, 0, 1]; [1, 1, 1, 0]; [1, 1, 1, 0]; [0, 0, 0, 1]];
test_left_edge_check = [[0, 0, 0, 0]; [1, 1, 1, 0]; [1, 1, 1, 0]; [0, 0, 0, 0]];
tested_left_edge = main_unit(test_left_edge);
assert(isequal(tested_left_edge, test_left_edge_check));

test_right_edge = [[1, 0, 0, 0]; [0, 1, 1, 1]; [0, 1, 1, 1]; [1, 0, 0, 0]];
test_right_edge_check = [[0, 0, 0, 0]; [0, 1, 1, 1]; [0, 1, 1, 1]; [0, 0, 0, 0]];
tested_right_edge = main_unit(test_right_edge);
assert(isequal(tested_right_edge, test_right_edge_check));





