function test_fatcurve
%TEST_FATCURVE Test the fatcurve function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help fatcurve
N=5;
phi = linspace(0, pi/5, N);
theta = zeros(1, N);
s = [theta; phi];
c = polar2cart(s)
r = 0.1
[X,Y,Z] = fatcurve(c,r)
