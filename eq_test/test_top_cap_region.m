function test_top_cap_region
%TEST_TOP_CAP_REGION Test the top_cap_region function

% Copyright 2025 Paul Leopardi.
% $Revision $ 1.12.2 $ $Date 2025-09-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help top_cap_region
region = top_cap_region(1, pi/6)
region = top_cap_region(2, pi/6)
region = top_cap_region(3, pi/6)
