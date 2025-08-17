function points = centres_of_regions(regions)
%CENTRES_OF_REGIONS Centre points of given regions
%
%Syntax
% points = centres_of_regions(regions);
%
%Description
% POINTS = CENTRES_OF_REGIONS(REGIONS) sets POINTS to be the centres of the
% regions in REGIONS, given as pairs of points in spherical polar coordinates.
% If REGIONS is a dim by 2 by N array, then the result POINTS is a dim by N
% array, in spherical polar coordinates.
%
%Examples
%
% >> dim = 2; N = 4;
% >> regions = eq_regions(dim, N);
% >> points = centres_of_regions(regions)
%
% points =
%
%          0    1.5708    4.7124         0
%          0    1.5708    1.5708    3.1416
%
% >> dim = 2; N = 10;
% >> regions = eq_regions(dim, N);
% >> points = centres_of_regions(regions)
%
% points =
%
%   Columns 1 through 7
%
%          0    0.7854    2.3562    3.9270    5.4978    1.5708    3.1416
%          0    1.1071    1.1071    1.1071    1.1071    2.0344    2.0344
%
%   Columns 8 through 10
%
%     4.7124         0         0
%     2.0344    2.0344    3.1416
%
% >> dim = 3; N = 6;
% >> regions = eq_regions(dim, N);
% >> points = centres_of_regions(regions)
%
% points =
%
%          0         0    1.5708    4.7124         0         0
%          0         0    1.5708    1.5708    3.1416         0
%          0    1.5708    1.5708    1.5708    1.5708    3.1416

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

tol = eps*2^5;
dim = size(regions,1);
N =   size(regions,3);
points = zeros(dim,N);
top = regions(:,1,:);
bot = regions(:,2,:);
zero_bot = abs(bot(1,:)) < tol;
bot(1,zero_bot) = 2*pi;
equal_bot = abs(bot(1,:) - top(1,:)) < tol;
bot(1,equal_bot) = top(1,equal_bot) + 2*pi;
twopi_bot = abs(bot(1,:) - top(1,:) - 2*pi) < tol;
points(1,twopi_bot) = 0;
points(1,~twopi_bot) = mod((bot(1,~twopi_bot) + top(1,~twopi_bot))/2, 2*pi);
for k = 2:dim
   pi_bot =   abs(bot(k,:) - pi) < tol;
   points(k,pi_bot) = pi;
   zero_top = abs(top(k,:)) < tol;
   points(k,zero_top) = 0;
   all_else = ~(zero_top | pi_bot);
   points(k,all_else) = mod((top(k,all_else) + bot(k,all_else))/2, pi);
end
%
% end function
