function test_euclidean_dist
%TEST_EUCLIDEAN_DIST Test the euclidean_dist function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help euclidean_dist
x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
y = [[0 -0.5 0.866]' [0 0.866 0.5]' [0 -0.866 -0.5]' [0 0.5 -0.866]']
distance = euclidean_dist(x,y)
