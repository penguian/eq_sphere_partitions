function test_eq_count_points_by_s2_region
%TEST_EQ_COUNT_POINTS_BY_S2_REGION Test the eq_count_points_by_s2_region function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_count_points_by_s2_region
points_s = eq_point_set_polar(2,8)
count_v = eq_count_points_by_s2_region(points_s, 8)
count_v = eq_count_points_by_s2_region(points_s, 5)
sum(count_v)
points_s = eq_point_set_polar(2,128);
count_v = eq_count_points_by_s2_region(points_s, 8)
sum(count_v)
count_v = eq_count_points_by_s2_region(points_s, 5)
sum(count_v)
