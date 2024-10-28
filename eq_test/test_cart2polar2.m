function test_cart2polar2
%TEST_CART2POLAR2 Test the cart2polar2 function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help cart2polar2
x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
s = cart2polar2(x)
