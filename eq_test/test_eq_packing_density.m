function test_eq_packing_density
%TEST_EQ_PACKING_DENSITY Test the eq_packing_density function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_packing_density
density = eq_packing_density(2,10)
density = eq_packing_density(3,1:6)
