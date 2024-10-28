function test_point_set_packing_density
%TEST_POINT_SET_PACKING_DENSITY Test the point_set_packing_density function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help point_set_packing_density
x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
density = point_set_packing_density(x)
