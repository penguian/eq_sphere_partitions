function c_caps = cap_colats(dim,N,c_polar,n_regions)
%CAP_COLATS Colatitudes of spherical caps enclosing cumulative sum of regions
%
%Syntax
% c_caps = cap_colats(dim,N,c_polar,n_regions);
%
%Description
% C_CAPS = CAP_COLATS(DIM,N,C_POLAR,N_REGIONS) determine C_CAPS, an increasing
% array of colatitudes of spherical caps which enclose the same area as that
% given by the cumulative sum of regions, given DIM, N, C_POLAR and N_REGIONS,
% where DIM is the dimension of the sphere, N is the number of regions in the
% equal area partition, C_POLAR is the colatitude of the North polar cap, and
% N_REGIONS is a list of the natural number of regions in each collar and the
% polar caps.
%
%Notes
% The length of both N_REGIONS and C_CAPS is N_COLLARS+2, where N_COLLARS is
% the number of collars. C_CAPS[1] is C_POLAR.
% C_CAPS[N_COLLARS+1] is Pi-C_POLAR. C_CAPS[N_COLLARS+2] is Pi.
%
%Examples
%
% >> dim = 2; N = 4;
% >> c_polar = polar_colat(dim, N); [s_cap, n_regions] = eq_caps(dim, N);
% >> c_caps = cap_colats(dim, N, c_polar, n_regions)
% >> assert(all(c_caps == s_cap));
%
% c_caps =
%
%     1.0472    2.0944    3.1416
%
% >> dim = 2; N = 10;
% >> c_polar = polar_colat(dim, N); [s_cap, n_regions] = eq_caps(dim, N);
% >> c_caps = cap_colats(dim, N, c_polar, n_regions)
% >> assert(all(c_caps == s_cap));
%
% c_caps =
%
%     0.6435    1.5708    2.4981    3.1416
%
% >> dim = 3; N = 6;
% >> c_polar = polar_colat(dim, N); [s_cap, n_regions] = eq_caps(dim, N);
% >> c_caps = cap_colats(dim, N, c_polar, n_regions)
% >> assert(all(c_caps == s_cap));
%
% c_caps =
%
%     0.9845    2.1571    3.1416

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2025-08-16 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

c_caps = zeros(size(n_regions));
c_caps(1) = c_polar;
ideal_region_area = area_of_ideal_region(dim,N);
n_collars = size(n_regions,2)-2;
subtotal_n_regions = 1;
for collar_n = 1:n_collars
    subtotal_n_regions = subtotal_n_regions+n_regions(1+collar_n);
    c_caps(collar_n+1) =sradius_of_cap(dim,subtotal_n_regions*ideal_region_area);
end
c_caps(1+n_collars+1) = pi;
%
% end function
