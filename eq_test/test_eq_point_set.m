function test_eq_point_set
%TEST_EQ_POINT_SET Test the eq_point_set function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_point_set
points_x = eq_point_set(2,4)
size(points_x)
