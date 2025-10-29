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
function surf_jet
%SURF_JET Set up extreme color values using COLORMAP JET
%
%Syntax
% surf_jet
%
%Examples
%
% >> surf_jet
% >> show_s2_sphere
%
%See also
% SHOW_S2_SPHERE, SHOW_S2_PARTITION

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

[X,Y,Z] = sphere(2);
surf(X*0,Y*0,Z*0,-ones(size(Z)),'FaceAlpha',0)
hold on
surf(X*0,Y*0,Z*0,ones(size(Z)),'FaceAlpha',0)
hold on
colormap jet
axis equal
grid off
axis off
if ~haslight(gca)
    camlight right
end
%
%end function
function result = x2eqarea(x)
%X2EQUAREA Equal area projection of 3D Euclidean points
%
%Syntax
% result = x2eqarea(x);
%
%Examples
%
% >> x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
%
% x =
%
%      0     0     0     0
%      0     1    -1     0
%      1     0     0    -1
%
% >> x2eqarea(points_x)
%
% ans =
%
%          0         0         0
%     1.4142   -1.4142         0
%
%See also
% X2STEREO, PROJECT_POINT_SET

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

dim = size(x,1)-1;
theta = acos(-x(dim+1,:));
r = (area_of_cap(dim,theta)/volume_of_ball(dim)).^(1/dim);
mask = (sin(theta) == 0);
scale = ones(dim,1)*(r(~mask)./sin(theta(~mask)));
eqarea = zeros(dim,size(x,2));
eqarea(:,~mask) = x(1:dim,~mask).*scale;
mask1 = (x(dim+1,:) == 1);
result = eqarea(:,~mask1);
% end function
function result = x2stereo(x)
%X2STEREO Stereographic projection of Euclidean points
%
%Syntax
% result = x2stereo(x);
%
%Examples
%
% >> x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
%
% x =
%
%      0     0     0     0
%      0     1    -1     0
%      1     0     0    -1
%
% >> x2stereo(points_x)
%
% ans =
%
%      0     0     0
%      1    -1     0
%
%See also
% X2EQAREA, PROJECT_POINT_SET

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

dim = size(x,1)-1;
mask = (x(dim+1,:) == 1);
scale = ones(dim,1)*(1-x(dim+1,~mask));
result = x(1:dim,~mask)./scale;
% end function
