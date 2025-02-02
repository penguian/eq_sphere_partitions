function test_point_set_energy_dist
%TEST_POINT_SET_ENERGY_DIST Test the point_set_energy_dist function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help point_set_energy_dist
x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
[energy,min_dist] = point_set_energy_dist(x)
