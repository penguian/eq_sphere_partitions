function test_eq_find_s2_region
%TEST_EQ_FIND_S2_REGION Test the eq_find_s2_region function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_find_s2_region
points_s = eq_point_set_polar(2,8)
r_idx = eq_find_s2_region(points_s,8)
r_idx = eq_find_s2_region(points_s,5)
