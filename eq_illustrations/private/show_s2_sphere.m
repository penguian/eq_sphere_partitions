function show_s2_sphere
%SHOW_S2_SPHERE Illustrate the unit sphere S^2
%
%Syntax
% show_s2_sphere
%
%Examples
%
% >> surf_jet
% >> show_s2_sphere
%
%See also
% SURF_JET, SHOW_S2_PARTITION

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

[X,Y,Z] = sphere(90);
surf(X,Y,Z,zeros(size(Z)),...
     'FaceColor','interp','FaceLighting','phong',...
     'EdgeColor','none','SpecularStrength',0)
axis equal
grid off
axis off
%
%end function
