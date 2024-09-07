function r_idx = eq_find_s2_region(s_point,N)
%EQ_FIND_S2_REGION Given a sequence of points, partition S^2 into N regions
%and find the index of the region containing each point.
%
%Syntax
% r_idx = eq_find_s2_region(s_point,N);
%
%Description
% r_idx = eq_find_s2_region(s_point,N) does the following:
% 1) partitions the unit sphere S^2 into a sequence of N regions of
%    equal area and small diameter;
% 2) for each point in s_point, determines which region contains the point;
% 3) sets r_idx to be an array of length size(s_point,2) containing the index of
%    the region corresponding to each point.
%
%Arguments
% s_point  Sequence of points on S^2, as a 2 x n_points array in spherical polar coordinates,
%          with longitude 0 <= s(1,p_idx) <= 2*pi, colatitude 0 <= s(2,p_idx) <= pi.
% N        Required number of regions, a positive integer.
%
%Examples
% > points_s = eq_point_set_polar(2,8)
% points_s =
%         0    0.5236    1.5708    2.6180    3.6652    4.7124    5.7596         0
%         0    1.5708    1.5708    1.5708    1.5708    1.5708    1.5708    3.1416
%
% > r_idx = eq_find_s2_region(points_s,5)
% r_idx =
%     1     2     2     3     3     4     4     5
%
%See also
% EQ_COUNT_POINTS_BY_S2_REGION, LOOKUP_S2_REGION

% Copyright 2024 Paul Leopardi
% $Revision 1.12 $ $Date 2024-09-07 $
% Copyright 2012 Paul Leopardi
% $Revision 1.11 $ $Date 2012-01-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

% Partition the sphere into N regions.
s_regions = eq_regions(2, N);
% Partition the sphere into nested spherical caps, making
% s_cap a sequence of cap colatitudes, and
% n_regions being the number of regions in each collar between
% successive colatitudes.
[s_cap, n_regions] = eq_caps(2, N);
% Set c_regions to be the cumulative number of regions within each nested cap.
c_regions = cumsum(n_regions);
% Set r_idx to be the sequence of region indices into s_regions, such that
% for all i from 1 to size(s_point, 2), the point s_point(i)
% is contained in the region s_regions(r_idx(i)).
r_idx = lookup_s2_region(s_point, s_regions, s_cap, c_regions);
