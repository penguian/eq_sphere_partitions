function test_ideal_region_list
%TEST_IDEAL_REGION_LIST Test the ideal_region_list function

% Copyright 2025 Paul Leopardi.
% $Revision $ 1.12.2 $ $Date 2025-08-17 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help ideal_region_list
dim = 2; N = 4;
c_polar = polar_colat(dim, N);
n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
r_regions = ideal_region_list(dim,N,c_polar,n_collars)
dim = 2; N = 10;
c_polar = polar_colat(dim, N);
n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
r_regions = ideal_region_list(dim,N,c_polar,n_collars)
dim = 2; N = 6;
c_polar = polar_colat(dim, N);
n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
r_regions = ideal_region_list(dim,N,c_polar,n_collars)
