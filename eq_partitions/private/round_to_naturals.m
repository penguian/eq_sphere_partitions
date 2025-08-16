function n_regions = round_to_naturals(N,r_regions)
%ROUND_TO_NATURALS Round off a given list of numbers of regions
%
%Syntax
% n_regions = round_to_naturals(N,r_regions);
%
%Description
% N_REGIONS = ROUND_TO_NATURALS(N,R_REGIONS) determines N_REGIONS,
% a list of the natural number of regions in each collar and the polar caps,
% given N, the number of regions in a sphere partion, and R_REGIONS, an
% array representing the ideal number of regions in each collar and the caps.
% The result N_REGIONS is as close as possible to R_REGIONS, using rounding.
%
%Notes
% The sum of both R_REGIONS and N_REGIONS is N. The length of both R_REGIONS
% and N_REGIONS is N_COLLARS+2, where N_COLLARS is the number of collars.
% N_REGIONS[1] is 1. N_REGIONS[N_COLLARS+2] is 1.
%
%Examples
%
% >> dim = 2; N = 4;
% >> c_polar = polar_colat(dim,N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars);
% >> n_regions = round_to_naturals(N,r_regions)
%
% n_regions =
%
%      1     2     1
%
% >> dim = 2; N = 10;
% >> c_polar = polar_colat(dim,N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars);
% >> n_regions = round_to_naturals(N,r_regions)
%
% n_regions =
%
%      1     4     4     1
%
% >> dim = 3; N = 6;
% >> c_polar = polar_colat(dim,N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars);
% >> n_regions = round_to_naturals(N,r_regions)
%
% n_regions =
%
%      1     4     1

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2025-08-16 $
% Added examples
% Copyright 2025 Paul Leopardi.
% $Revision 1.12.1 $ $Date 2025-08-09 $
% Added an assertion.
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-05-05 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_regions = r_regions;
discrepancy = 0;
for zone_n = 1:size(r_regions,2)
    n_regions(zone_n) = round(r_regions(zone_n)+discrepancy);
    discrepancy = discrepancy+r_regions(zone_n)-n_regions(zone_n);
end
assert(sum(n_regions)==N,'Sum of result n_regions does not equal N==%g',N)
%
% end function
