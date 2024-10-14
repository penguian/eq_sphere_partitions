function test_point_set_dist_coeff
%TEST_POINT_SET_DIST_COEFF Test the point_set_dist_coeff function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help point_set_dist_coeff
x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
coeff = point_set_dist_coeff(x)
