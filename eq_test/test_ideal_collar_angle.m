function test_ideal_collar_angle
%TEST_IDEAL_COLLAR_ANGLE Test the ideal_collar_angle function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help ideal_collar_angle
angle = ideal_collar_angle(2,10)
angle = ideal_collar_angle(3,1:6)
