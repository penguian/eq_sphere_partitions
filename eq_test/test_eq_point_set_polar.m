function test_eq_point_set_polar
%TEST_EQ_POINT_SET_POLAR Test the eq_point_set_polar function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_point_set_polar
points_s = eq_point_set_polar(2,4)
size(points_s)
