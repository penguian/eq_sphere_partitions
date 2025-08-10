function region = bot_cap_region(dim,a_cap)
%BOT_CAP_REGION South polar (bottom) cap region of EQ partition
%
% An array of two points representing the bottom cap of radius a_cap as a region.
%
% region = bot_cap_region(dim,a_cap);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

if dim == 1
    region = [2*pi-a_cap,2*pi];
else
    sphere_region_1 = sphere_region(dim-1);
    region = [[sphere_region_1(:,1); pi-a_cap],[sphere_region_1(:,2); pi]];
end
%
% end function
function c_caps = cap_colats(dim,N,c_polar,n_regions)
%CAP_COLATS Colatitudes of spherical caps enclosing cumulative sum of regions
%
% Given dim, N, c_polar and n_regions, determine c_caps,
% an increasing list of colatitudes of spherical caps which enclose the same area
% as that given by the cumulative sum of regions.
% The number of elements is n_collars+2.
% c_caps[1] is c_polar.
% c_caps[n_collars+1] is Pi-c_polar.
% c_caps[n_collars+2] is Pi.
%
% c_caps = cap_colats(dim,N,c_polar,n_regions);

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

% end function
function points = centres_of_regions(regions)
%CENTRES_OF_REGIONS Centre points of given regions
%
% points = centres_of_regions(regions);

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
% end function
function offset = circle_offset(n_top,n_bot,extra_twist)
%CIRCLE_OFFSET Try to maximize minimum distance of center points for S^2 collars
%
% Given n_top and n_bot, calculate an offset.
%
% The values n_top and n_bot represent the numbers of
% equally spaced points on two overlapping circles.
% The offset is given in multiples of whole rotations, and
% consists of three parts;
% 1) Half the difference between a twist of one sector on each of bottom and top.
% This brings the centre points into alignment.
% 2) A rotation which will maximize the minimum angle between
% points on the two circles.
% 3) An optional extra twist by a whole number of sectors on the second circle.
% The extra twist is added so that the location of
% the minimum angle  between circles will
% progressively twist around the sphere with each collar.
%
% offset = circle_offset(n_top,n_bot,extra_twist);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

if nargin < 3
    extra_twist = false;
end
offset = (1/n_bot - 1/n_top)/2 + gcd(n_top,n_bot)/(2*n_top*n_bot);
if extra_twist
    twist = 6;
    offset = offset + twist/n_bot;
end
% end function
function r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%IDEAL_REGION_LIST The ideal real number of regions in each zone
%
% List the ideal real number of regions in each collar, plus the polar caps.
%
% Given dim, N, c_polar and n_collars, determine r_regions,
% a list of the ideal real number of regions in each collar,
% plus the polar caps.
% The number of elements is n_collars+2.
% r_regions[1] is 1.
% r_regions[n_collars+2] is 1.
% The sum of r_regions is N.
%
% r_regions = ideal_region_list(dim,N,c_polar,n_collars);

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

% end function
function n_collars = num_collars(N,c_polar,a_ideal)
%NUM_COLLARS The number of collars between the polar caps
%
% Given N, an ideal angle, and c_polar,
% determine n_collars, the number of collars between the polar caps.
%
%  n_collars = num_collars(N,c_polar,a_ideal);

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_collars = zeros(size(N));
enough = (N > 2) & (a_ideal > 0);
n_collars(enough) = max(1,round((pi-2*c_polar(enough))./a_ideal(enough)));
%
% end function
function c_polar = polar_colat(dim, N)
%POLAR_COLAT The colatitude of the North polar (top) spherical cap
%
% Given dim and N, determine the colatitude of the North polar spherical cap.
%
% c_polar = polar_colat(dim, N);

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
function R = rot3(axis, angle)
%ROT3 R^3 rotation about a coordinate axis
%
%Syntax
% R = rot3(axis, angle);
%
%Description
% R = ROT3(AXIS, ANGLE) sets R to be the 3 by 3 rotation matrix corresponding to
% AXIS and ANGLE. Use this to create rotation matrices from Euler angles.
%
% AXIS must be 1, 2, or 3.
% ANGLE must be a real number.
%
%Examples
% > r=rot3(1,pi/6)
% r =
%     1.0000         0         0
%          0    0.8660   -0.5000
%          0    0.5000    0.8660
%
% > r=rot3(2,pi/6)
% r =
%     0.8660         0   -0.5000
%          0    1.0000         0
%     0.5000         0    0.8660
%
% > r=rot3(3,pi/6)
% r =
%     0.8660   -0.5000         0
%     0.5000    0.8660         0
%          0         0    1.0000

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

c = cos(angle);
s = sin(angle);
switch axis
case 1
    R = [...
            1, 0, 0; ...
            0, c, -s; ...
            0, s, c ];
case 2
    R = [...
            c, 0, -s; ...
            0, 1, 0; ...
            s, 0, c];
case 3
    R = [...
            c, -s, 0; ...
            s, c,  0; ...
            0, 0,  1];
end
%
%end function

function n_regions = round_to_naturals(N,r_regions)
%ROUND_TO_NATURALS Round off a given list of numbers of regions
%
% Given N and r_regions, determine n_regions,
% a list of the natural number of regions in each collar and the polar caps.
% This list is as close as possible to r_regions, using rounding.
% The number of elements is n_collars+2.
% n_regions[1] is 1.
% n_regions[n_collars+2] is 1.
% The sum of n_regions is N.
%
% n_regions = round_to_naturals(N,r_regions);

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.1 $ $Date 2025-08-09 $
% Add an assertion.
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-05-05 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_regions = r_regions;
discrepancy = 0;
for zone_n = 1:size(r_regions,2)
    n_regions(zone_n) = round(r_regions(zone_n)+discrepancy);
    discrepancy = discrepancy+r_regions(zone_n)-n_regions(zone_n);
end
assert(sum(n_regions)==N,'Sum of result n_regions does not equal N==%g',N)
%
% end function
function rotation = s2_offset(points_1)
%S2_OFFSET Experimental offset rotation of S^2
%
%Syntax
% rotation = s2_offset(points_1);
%
%Description
% ROTATION = S2_OFFSET(POINTS_1) sets ROTATION to be an R^3 rotation matrix which
% rotates the north pole of S^2 to a point specified by the points of POINTS_1.
%
% POINTS_1 must be a 2 by M matrix, representing M points of S^2 in spherical
% polar coordinates, with M a positive integer.
%
%Examples
% > s
% s =
%          0    0.7854    2.3562    3.9270    5.4978         0
%          0    1.5708    1.5708    1.5708    1.5708    3.1416
%
% > r=s2_offset(s)
% r =
%     0.0000    0.7071    0.7071
%    -1.0000    0.0000         0
%    -0.0000   -0.7071    0.7071
%
%See also
% ROT3, CIRCLE_OFFSET

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_in_collar = size(points_1,2);
if n_in_collar > 2
    if (n_in_collar > 3) && (points_1(2,2) == points_1(2,3))
        a_3 = (points_1(1,2) + points_1(1,3))/2;
    else
        a_3 = points_1(1,2) + pi;
    end
    a_2 = points_1(2,2)/2;
else
    a_3 = 0;
    a_2 = pi/2;
end
rotation = rot3(2,-a_2)*rot3(3,-a_3);
%
%end function
function region = sphere_region(dim)
%SPHERE_REGION The sphere represented as a single region of an EQ partition
%
% An array of two points representing S^dim as a region.
%
% region = sphere_region(dim);

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
function region = top_cap_region(dim,a_cap)
%TOP_CAP_REGION North polar (top) cap region of EQ partition
%
% An array of two points representing the top cap of radius a_cap as a region.
%
% region = top_cap_region(dim,a_cap);

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
% end function
