function test_eq_area_error
%TEST_EQ_AREA_ERROR Test the eq_area_error function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help eq_area_error
[total_error, max_error] = eq_area_error(2,10)
[total_error, max_error] = eq_area_error(3,1:6)
