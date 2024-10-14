function test_calc_dist_coeff
%TEST_CALC_DIST_COEFF Test the calc_dist_coeff function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help calc_dist_coeff
N = 2:6
dist = eq_min_dist(2,N)
calc_dist_coeff(2,N,dist)
