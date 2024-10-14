function test_eq_point_set_property
%TEST_EQ_POINT_SET_PROPERTY Test the eq_point_set_property function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_point_set_property
dist = eq_point_set_property(@point_set_min_dist,2,10)
