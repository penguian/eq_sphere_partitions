function test_eq_count_points_by_s2_region
%TEST_EQ_COUNT_POINTS_BY_S2_REGION Test the eq_count_points_by_s2_region function

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.3 $ $Date 2026-04-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

points_s = eq_point_set_polar(2, 8);

% Test N=8
count_v = eq_count_points_by_s2_region(points_s, 8);
assert(isequal(count_v, ones(1, 8)), 'N=8 count mismatch');

% Test N=5
count_v = eq_count_points_by_s2_region(points_s, 5);
assert(isequal(count_v, [1, 2, 2, 2, 1]), 'N=5 count mismatch');
assert(sum(count_v) == 8, 'N=5 sum mismatch');

points_s = eq_point_set_polar(2, 128);

% Test N=8 with larger set
count_v = eq_count_points_by_s2_region(points_s, 8);
assert(isequal(count_v, [19, 15, 14, 17, 15, 14, 15, 19]), 'N=128, N_regions=8 count mismatch');
assert(sum(count_v) == 128, 'N=128, N_regions=8 sum mismatch');

% Test N=5 with larger set
count_v = eq_count_points_by_s2_region(points_s, 5);
assert(isequal(count_v, [19, 29, 32, 29, 19]), 'N=128, N_regions=5 count mismatch');
assert(sum(count_v) == 128, 'N=128, N_regions=5 sum mismatch');

fprintf('test_eq_count_points_by_s2_region: PASS\n');
