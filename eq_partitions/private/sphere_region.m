function region = sphere_region(dim)
%SPHERE_REGION The sphere represented as a single region of an EQ partition
%
%Syntax
% region = sphere_region(dim);
%
%Description
% REGION = SPHERE_REGION(DIM) Sets REGION to be array of two points in
% spherical polar coordinates representing the sphere S^DIM as a single region.
%
%Examples
%
% >> region = sphere_region(1)
%
% region =
%
%          0    6.2832
%
% >> region = sphere_region(2)
%
% region =
%
%          0    6.2832
%          0    3.1416
%
% >> region = sphere_region(3)
%
% region =
%
%          0    6.2832
%          0    3.1416
%          0    3.1416
%
%See also
% BOT_CAP_REGION, TOP_CAP_REGION

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
    region = [0,2*pi];
else
    sphere_region_1 = sphere_region(dim-1);
    region = [[sphere_region_1(:,1); 0],[sphere_region_1(:,2); pi]] ;
end
%
% end function
