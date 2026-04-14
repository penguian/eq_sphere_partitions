function test_in_s2_region
%TEST_IN_S2_REGION Test the in_s2_region function

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.3 $ $Date 2026-04-14 $
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

points_s = eq_point_set_polar(2, 8);
s_regions = eq_regions(2, 5);

% Check region 3 (index 3 in MATLAB)
region = s_regions(:, :, 3);
result = in_s2_region(points_s, region);

% Expected: points 4 and 5 (from eq_find_s2_region logic)
expected = [false, false, false, true, true, false, false, false];
assert(isequal(result, expected), 'in_s2_region mismatch for region 3');

fprintf('test_in_s2_region: PASS\n');
