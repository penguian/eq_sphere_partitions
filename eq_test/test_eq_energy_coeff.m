function test_eq_energy_coeff
%TEST_EQ_ENERGY_COEFF Test the eq_energy_coeff function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_energy_coeff
coeff = eq_energy_coeff(2,10)
coeff = eq_energy_coeff(3,1:6)
coeff = eq_energy_coeff(2,1:6,0)
