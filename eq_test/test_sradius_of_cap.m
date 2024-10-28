function test_sradius_of_cap
%TEST_SRADIUS_OF_CAP Test the sradius_of_cap function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help sradius_of_cap
s_cap = sradius_of_cap(2,area_of_sphere(2)/2)
s_cap = sradius_of_cap(3,(0:4)*area_of_sphere(3)/4)
