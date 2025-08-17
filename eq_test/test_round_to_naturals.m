function test_round_to_naturals
%TEST_ROUND_TO_NATURALS Test the round_to_naturals function

% Copyright 2025 Paul Leopardi.
% $Revision $ 1.12.2 $ $Date 2025-08-17 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help round_to_naturals
dim = 2; N = 4;
c_polar = polar_colat(dim,N);
n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
r_regions = ideal_region_list(dim,N,c_polar,n_collars);
n_regions = round_to_naturals(N,r_regions)
dim = 2; N = 10;
c_polar = polar_colat(dim,N);
n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
r_regions = ideal_region_list(dim,N,c_polar,n_collars);
n_regions = round_to_naturals(N,r_regions)
dim = 3; N = 6;
c_polar = polar_colat(dim,N);
n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
r_regions = ideal_region_list(dim,N,c_polar,n_collars);
n_regions = round_to_naturals(N,r_regions)
