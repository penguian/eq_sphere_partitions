function test_num_collars
%TEST_NUM_COLLARS Test the num_collars function

% Copyright 2025 Paul Leopardi.
% $Revision $ 1.12.2 $ $Date 2025-09-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help num_collars
dim = 2; N = 4;
c_polar = polar_colat(dim, N);
a_ideal = ideal_collar_angle(dim, N);
n_collars = num_collars(N,c_polar,a_ideal)
dim = 2; N = 10;
c_polar = polar_colat(dim, N);
a_ideal = ideal_collar_angle(dim, N);
n_collars = num_collars(N,c_polar,a_ideal)
dim = 3; N = 6;
c_polar = polar_colat(dim, N);
a_ideal = ideal_collar_angle(dim, N);
n_collars = num_collars(N,c_polar,a_ideal)
