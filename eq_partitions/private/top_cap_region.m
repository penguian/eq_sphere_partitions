function region = top_cap_region(dim,a_cap)
%TOP_CAP_REGION North polar (top) cap region of EQ partition
%
%Syntax
% region = top_cap_region(dim,a_cap);
%
%Description
% REGION = TOP_CAP_REGION(DIM,A_CAP) Sets REGION to be an array of two
% points in spherical polar coordinates representing the North polar
% (top) cap of spherical radius A_CAP as a region of the sphere S^DIM.
%
%Examples
%
% >> region = top_cap_region(1, pi/6)

% region =
%
%          0    0.5236
%
% >> region = top_cap_region(2, pi/6)
%
% region =
%
%          0    6.2832
%          0    0.5236
%
% >> region = top_cap_region(3, pi/6)
%
% region =
%
%          0    6.2832
%          0    3.1416
%          0    0.5236
%
%See also
% BOT_CAP_REGION, SPHERE_REGION

% Copyright 2025 Paul Leopardi
% $Revision 1.12.2 $ $Date 2025-08-16 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

if dim == 1
    region = [0,a_cap];
else
    sphere_region_1 = sphere_region(dim-1);
    region = [[sphere_region_1(:,1); 0], [sphere_region_1(:,2); a_cap]];
end
%
% end function
