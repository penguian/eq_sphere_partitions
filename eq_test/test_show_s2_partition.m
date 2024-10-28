function test_show_s2_partition
%TEST_SHOW_S2_PARTITION Test the show_s2_partition function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help show_s2_partition
show_s2_partition(10)
frames = show_s2_partition(9,'offset','extra')
show_s2_partition(99,'points','hide')
