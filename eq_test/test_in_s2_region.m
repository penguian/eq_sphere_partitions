function test_in_s2_region
%TEST_IN_S2_REGION Test the in_s2_region function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help in_s2_region
points_s = eq_point_set_polar(2,8)
s_regions = eq_regions(2,5);
region = s_regions(:,:,3)
result = in_s2_region(points_s, region)
