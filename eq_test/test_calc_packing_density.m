function test_calc_packing_density
%TEST_CALC_PACKING_DENSITY Test the calc_packing_density function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help calc_packing_density
N = 2:6
dist = eq_min_dist(2,N)
density = calc_packing_density(2,N,dist)
