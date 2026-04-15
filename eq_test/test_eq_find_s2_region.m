function test_eq_find_s2_region
%TEST_EQ_FIND_S2_REGION Test the eq_find_s2_region function

% Copyright 2024-2026 Paul Leopardi.
% $Revision 1.12.3 $ $Date 2026-04-14 $
% $Revision 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

points_s = eq_point_set_polar(2, 8);

% Test N=8
r_idx = eq_find_s2_region(points_s, 8);
assert(isequal(r_idx, 1:8), 'N=8 find mismatch');

% Test N=5
r_idx = eq_find_s2_region(points_s, 5);
assert(isequal(r_idx, [1, 2, 2, 3, 3, 4, 4, 5]), 'N=5 find mismatch');

fprintf('test_eq_find_s2_region: PASS\n');
