function test_calc_energy_coeff
%TEST_CALC_ENERGY_COEFF Test the calc_energy_coeff function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help calc_energy_coeff
N = 2:6
energy = eq_energy_dist(2,N,0)
coeff = calc_energy_coeff(2,N,0,energy)
