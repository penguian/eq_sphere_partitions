function test_eq_diam_coeff
%TEST_EQ_DIAM_COEFF Test the eq_diam_coeff function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_diam_coeff
bound_coeff = eq_diam_coeff(2,10)
[bound_coeff,vertex_coeff]=eq_diam_coeff(3,1:6)
