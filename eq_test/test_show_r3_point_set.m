function test_show_r3_point_set
%TEST_SHOW_R3_POINT_SET Test the show_r3_point_set function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help show_r3_point_set
points_x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
show_r3_point_set(points_x,'sphere','hide')
show_r3_point_set(points_x,'sphere','show')
