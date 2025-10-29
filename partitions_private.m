function region = bot_cap_region(dim,a_cap)
%BOT_CAP_REGION South polar (bottom) cap region of EQ partition
%
%Syntax
% region = bot_cap_region(dim,a_cap);
%
%Description
% REGION = BOT_CAP_REGION(DIM,A_CAP) Sets REGION to be an array of two
% points in spherical polar coordinates representing the South polar
% (bottom) cap of spherical radius A_CAP as a region of the sphere
% S^DIM.
%
%Examples
%
% >> region = bot_cap_region(1, pi/6)
%
% region =
%
%     5.7596    6.2832
%
% >> region = bot_cap_region(2, pi/6)
%
% region =
%
%          0    6.2832
%     2.6180    3.1416
%
% >> region = bot_cap_region(3, pi/6)
%
% region =
%
%          0    6.2832
%          0    3.1416
%     2.6180    3.1416
%
%See also
% TOP_CAP_REGION, SPHERE_REGION

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
function offset = circle_offset(n_top,n_bot,extra_twist)
%CIRCLE_OFFSET Try to maximize minimum distance of center points for S^2 collars
%
%Syntax
% offset = circle_offset(n_top,n_bot,extra_twist);
%
%Description
% OFFSET = CIRCLE_OFFSET(N_TOP,N_BOT,EXTRA_TWIST) calculates the offset OFFSET
% that maximizes the minimum distance between two sets of points on a circle,
% determined by the numbers N_TOP and N_BOT.
%
%Notes
% The values N_TOP and N_BOT represent the numbers of points in each of two
% equally spaced sets on the unit circle. The offset is given in multiples of
% whole rotations, and consists of three parts:
% 1) Half the difference between a twist of one sector on each of bottom and top.
% This brings the centre points into alignment.
% 2) A rotation which will maximize the minimum angle between points on the two
% circles.
% 3) An optional extra twist by a whole number of sectors on the second circle.
% The extra twist is added so that the location of the minimum angle between
% circles will progressively twist around the sphere with each collar.
%
%Examples
%
% >> offset = circle_offset(3, 4)
%
% offset =
%
%    6.9389e-18
%
% >> offset = circle_offset(3, 4, 1)
%
% offset =
%
%     1.5000
%
% >> offset = circle_offset(7, 11)
%
% offset =
%
%    -0.0195

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

if nargin < 3
    extra_twist = false;
end
offset = (1/n_bot - 1/n_top)/2 + gcd(n_top,n_bot)/(2*n_top*n_bot);
if extra_twist
    twist = 6;
    offset = offset + twist/n_bot;
end
%
% end function
function r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%IDEAL_REGION_LIST The ideal real number of regions in each zone
%
%Syntax
% r_regions = ideal_region_list(dim,N,c_polar,n_collars);
%
%Description
% R_REGIONS = IDEAL_REGION_LIST(DIM,N,C_POLAR,N_COLLARS) determines R_REGIONS,
% an array of the ideal real number of regions in each collar plus the polar
% caps, given DIM, N, C_POLAR and N_COLLARS, where DIM is the dimension of the
% sphere S^DIM, N is the number of regions in the equal area partition, C_POLAR
% is the colatitude of the North polar (top) spherical cap, and N_COLLARS is
% the number of collars in the partition.
%
%Notes
% The length of R_REGIONS is N_COLLARS+2. R_REGIONS[1] is 1.
% R_REGIONS[N_COLLARS+2] is 1. The sum of R_REGIONS is N.
%
%Examples
% >> dim = 2; N = 4;
% >> c_polar = polar_colat(dim, N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%
% r_regions =
%
%     1.0000    2.0000    1.0000
%
% >> dim = 2; N = 10;
% >> c_polar = polar_colat(dim, N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%
% r_regions =
%
%     1.0000    4.0000    4.0000    1.0000
%
% >> dim = 2; N = 6;
% >> c_polar = polar_colat(dim, N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars)
%
% r_regions =
%
%     1.0000    4.0000    1.0000

%See also
% IDEAL_COLLAR_ANGLE, NUM_COLLARS, POLAR_COLAT

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
%
% end function
function n_collars = num_collars(N,c_polar,a_ideal)
%NUM_COLLARS The number of collars between the polar caps
%
%Syntax
%  n_collars = num_collars(N,c_polar,a_ideal);
%
%Description
% N_COLLARS = NUM_COLLARS(N,C_POLAR,A_IDEAL) determines N_COLLARS, the number
% of collars between the polar caps of an equal area partition of the sphere
% into N regions, given the North polar colatitude C_POLAR, and the ideal
% collar angle A_IDEAL.
%
%Examples
%
% >> dim = 2; N = 4;
% >> c_polar = polar_colat(dim, N);
% >> a_ideal = ideal_collar_angle(dim, N);
% >> n_collars = num_collars(N,c_polar,a_ideal)
%
% n_collars =
%
%      1
%
% >> dim = 2; N = 10;
% >> c_polar = polar_colat(dim, N);
% >> a_ideal = ideal_collar_angle(dim, N);
% >> n_collars = num_collars(N,c_polar,a_ideal)
%
% n_collars =
%
%      2
%
% >> dim = 3; N = 6;
% >> c_polar = polar_colat(dim, N);
% >> a_ideal = ideal_collar_angle(dim, N);
% >> n_collars = num_collars(N,c_polar,a_ideal)
%
% n_collars =
%
%      1
%
%See also
% IDEAL_COLLAR_ANGLE, IDEAL_REGION_LIST, POLAR_COLAT

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

n_collars = zeros(size(N));
enough = (N > 2) & (a_ideal > 0);
n_collars(enough) = max(1,round((pi-2*c_polar(enough))./a_ideal(enough)));
%
% end function
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
%
% >> r = rot3(1,pi/6)
%
% r =
%     1.0000         0         0
%          0    0.8660   -0.5000
%          0    0.5000    0.8660
%
% >> r = rot3(2, pi/6)
%
% r =
%     0.8660         0   -0.5000
%          0    1.0000         0
%     0.5000         0    0.8660
%
% >> r = rot3(3, pi/6)
%
% r =
%     0.8660   -0.5000         0
%     0.5000    0.8660         0
%          0         0    1.0000

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
%Syntax
% n_regions = round_to_naturals(N,r_regions);
%
%Description
% N_REGIONS = ROUND_TO_NATURALS(N,R_REGIONS) determines N_REGIONS,
% a list of the natural number of regions in each collar and the polar caps,
% given N, the number of regions in a sphere partition, and R_REGIONS, an
% array representing the ideal number of regions in each collar and the caps.
% The result N_REGIONS is as close as possible to R_REGIONS, using rounding.
%
%Notes
% The sum of both R_REGIONS and N_REGIONS is N. The length of both R_REGIONS
% and N_REGIONS is N_COLLARS+2, where N_COLLARS is the number of collars.
% N_REGIONS[1] is 1. N_REGIONS[N_COLLARS+2] is 1.
%
%Examples
%
% >> dim = 2; N = 4;
% >> c_polar = polar_colat(dim,N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars);
% >> n_regions = round_to_naturals(N,r_regions)
%
% n_regions =
%
%      1     2     1
%
% >> dim = 2; N = 10;
% >> c_polar = polar_colat(dim,N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars);
% >> n_regions = round_to_naturals(N,r_regions)
%
% n_regions =
%
%      1     4     4     1
%
% >> dim = 3; N = 6;
% >> c_polar = polar_colat(dim,N);
% >> n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
% >> r_regions = ideal_region_list(dim,N,c_polar,n_collars);
% >> n_regions = round_to_naturals(N,r_regions)
%
% n_regions =
%
%      1     4     1

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2025-08-16 $
% Added examples
% Copyright 2025 Paul Leopardi.
% $Revision 1.12.1 $ $Date 2025-08-09 $
% Added an assertion.
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
