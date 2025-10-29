function pseudo_region = pseudo_region_for_diam(region)
%PSEUDO_REGION_FOR_DIAM Two points which maximize the vertex diameter of a region
%
%Syntax
% pseudo_region = pseudo_region_for_diam(region)
%
%Examples
%
% >> regions = eq_regions(2,10);
% >> region = regions(:,:,2)
%
% region =
%
%          0    1.5708
%     0.6435    1.5708
%
% >> pseudo_region = pseudo_region_for_diam(region)
%
% pseudo_region =
%
%          0    1.5708
%     0.6435    1.5708
%
%See also
% EXPAND_REGION_FOR_DIAM, MAX_VERTEX_DIAM_OF_REGIONS, EQ_VERTEX_DIAM

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2025-10-29 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

tol = eps*2^5;
phi_top = region(1,1);
phi_bot = region(1,2);
if phi_bot == 0
    phi_bot = 2*pi;
end
if (mod(phi_bot - phi_top, 2*pi) < tol) || (mod(phi_bot - phi_top, 2*pi) > pi)
    phi_bot = phi_top + pi;
end
pseudo_region = region;
pseudo_region(1,2) = phi_bot;
%
% end function
