function test_cap_colats
%TEST_CAP_COLATS Test the cap_colats function

% Copyright 2025 Paul Leopardi.
% $Revision $ 1.12.2 $ $Date 2025-09-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help cap_colats
dim = 2; N = 4;
c_polar = polar_colat(dim, N); [s_cap, n_regions] = eq_caps(dim, N);
c_caps = cap_colats(dim, N, c_polar, n_regions)
assert(all(c_caps == s_cap));
dim = 2; N = 10;
c_polar = polar_colat(dim, N); [s_cap, n_regions] = eq_caps(dim, N);
c_caps = cap_colats(dim, N, c_polar, n_regions)
assert(all(c_caps == s_cap));
dim = 3; N = 6;
c_polar = polar_colat(dim, N); [s_cap, n_regions] = eq_caps(dim, N);
c_caps = cap_colats(dim, N, c_polar, n_regions)
assert(all(c_caps == s_cap));
