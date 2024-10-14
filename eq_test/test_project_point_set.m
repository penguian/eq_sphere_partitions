function test_project_point_set
%TEST_PROJECT_POINT_SET Test the project_point_set function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help project_point_set
x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
project_point_set(x)
