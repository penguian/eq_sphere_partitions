function test_sph2euc_dist
%TEST_SPH2EUC_DIST Test the sph2euc_dist function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help sph2euc_dist
e = sph2euc_dist(pi)
e = sph2euc_dist(0:pi/4:pi)
