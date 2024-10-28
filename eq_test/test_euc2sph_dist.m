function test_euc2sph_dist
%TEST_EUC2SPH_DIST Test the euc2sph_dist function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help euc2sph_dist
s = euc2sph_dist(2)
s = euc2sph_dist(0:0.5:2)
s = euc2sph_dist(-2)
