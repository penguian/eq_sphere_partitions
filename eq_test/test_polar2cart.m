function test_polar2cart
%TEST_POLAR2CART Test the polar2cart function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help polar2cart
s = [[0 0]' [pi/2 pi/2]' [3*pi/2 pi/2]' [0 pi]']
x = polar2cart(s)
