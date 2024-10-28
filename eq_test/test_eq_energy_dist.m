function test_eq_energy_dist
%TEST_EQ_ENERGY_DIST Test the eq_energy_dist function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_energy_dist
energy = eq_energy_dist(2,10)
[energy,dist] = eq_energy_dist(3,1:6,0)
[energy,dist] = eq_energy_dist(3,100,1,'offset','extra')
