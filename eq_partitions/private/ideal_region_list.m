function r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%IDEAL_REGION_LIST The ideal real number of regions in each zone
%
%Syntax
% r_regions = ideal_region_list(dim,N,c_polar,n_collars);
%
%Description
% R_REGIONS = IDEAL_REGION_LIST(DIM,N,C_POLAR,N_COLLARS) determines R_REGIONS,
% an array of the ideal real number of regions in each collar plus the polar
% caps, given DIM, N, C_POLAR and N_COLLARS, where DIM is the dimension of the
% sphere S^DIM, N is the number of regions in the equal area partition, C_POLAR
% is the colatitude of the North polar (top) spherical cap, and N_COLLARS is
% the number of collars in the partition.
%
%Notes
% The length of R_REGIONS is N_COLLARS+2. R_REGIONS[1] is 1.
% R_REGIONS[N_COLLARS+2] is 1. The sum of R_REGIONS is N.
%
%Examples
% >> dim = 2; N = 4;
% >> c_polar = polar_colat(dim, N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%
% r_regions =
%
%     1.0000    2.0000    1.0000
%
% >> dim = 2; N = 10;
% >> c_polar = polar_colat(dim, N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%
% r_regions =
%
%     1.0000    4.0000    4.0000    1.0000
%
% >> dim = 2; N = 6;
% >> c_polar = polar_colat(dim, N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%
% r_regions =
%
%     1.0000    4.0000    1.0000

%See also
% IDEAL_COLLAR_ANGLE, NUM_COLLARS, POLAR_COLAT

% Copyright 2025 Paul Leopardi
% $Revision 1.12.2 $ $Date 2025-08-17 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

r_regions = zeros(1,2+n_collars);
r_regions(1) = 1;
if n_collars > 0
    %
    % Based on n_collars and c_polar, determine a_fitting,
    % the collar angle such that n_collars collars fit between the polar caps.
    %
    a_fitting = (pi-2*c_polar)/n_collars;
    ideal_region_area = area_of_ideal_region(dim,N);
    for collar_n = 1:n_collars
        ideal_collar_area = area_of_collar(dim, c_polar+(collar_n-1)*a_fitting, c_polar+collar_n*a_fitting);
        r_regions(1+collar_n) = ideal_collar_area / ideal_region_area;
    end
end
r_regions(2+n_collars) = 1;
%
% end function
