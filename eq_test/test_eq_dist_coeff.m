function test_eq_dist_coeff
%TEST_EQ_DIST_COEFF Test the eq_dist_coeff function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_dist_coeff
coeff = eq_dist_coeff(2,10)
coeff = eq_dist_coeff(3,1:6)
