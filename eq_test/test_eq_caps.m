function test_eq_caps
%TEST_EQ_CAPS Test the eq_caps function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_caps
[s_cap,n_regions] = eq_caps(2,10)
[s_cap,n_regions] = eq_caps(3,6)
