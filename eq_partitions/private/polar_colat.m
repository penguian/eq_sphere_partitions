function c_polar = polar_colat(dim, N)
%POLAR_COLAT The colatitude of the North polar (top) spherical cap
%
%Syntax
% c_polar = polar_colat(dim, N);
%
%Description
% C_POLAR = POLAR_COLAT(DIM, N) determines the colatitude of the North polar
% spherical cap of the partition of the sphere S^DIM into N equal regions.
%
%Examples
%
% >> c_polar = polar_colat(2, 4)
%
% c_polar =
%
%     1.0472
%
% >> c_polar = polar_colat(2, 10)
%
% c_polar =
%
%     0.6435
%
% >> c_polar = polar_colat(3, 6)
%
% c_polar =
%
%     0.9845
%
%See also
% IDEAL_COLLAR_ANGLE, IDEAL_REGION_LIST, NUM_COLLARS

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

enough = N > 2;
c_polar(N == 1) = pi;
c_polar(N == 2) = pi/2;
c_polar(enough) = sradius_of_cap(dim,area_of_ideal_region(dim,N(enough)));
%
% end function
