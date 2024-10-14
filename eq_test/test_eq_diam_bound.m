function test_eq_diam_bound
%TEST_EQ_DIAM_BOUND Test the eq_diam_bound function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_diam_bound
diam_bound = eq_diam_bound(2,10)
diam_bound = eq_diam_bound(3,1:6)
