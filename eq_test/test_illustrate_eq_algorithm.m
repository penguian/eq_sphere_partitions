function test_illustrate_eq_algorithm
%TEST_ILLUSTRATE_EQ_ALGORITHM Test the illustrate_eq_algorithm function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help illustrate_eq_algorithm
illustrate_eq_algorithm(3,99)
illustrate_eq_algorithm(3,99,'offset','extra','proj','eqarea')
illustrate_eq_algorithm(3,99,'proj','eqarea','points','hide')
