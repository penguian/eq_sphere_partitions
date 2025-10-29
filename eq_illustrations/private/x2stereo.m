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
