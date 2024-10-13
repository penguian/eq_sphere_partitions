function count_v = eq_count_points_by_s2_region(s_point,N)
%EQ_COUNT_POINTS_BY_S2_REGION Given a set of points, count points in each of N regions of S^2
%
%Syntax
% count_v = eq_count_points_by_s2_region(s_point,N)
%
%Description
% count_v = eq_count_points_by_s2_region(s_point,N) does the following:
% 1) partitions the unit sphere S^2 into a sequence of N regions of
%    equal area and small diameter;
% 2) for each point in s_point, determines which region contains the point;
% 3) sets count_v to be an array of length N containing the number of points of
%    s_point contained in each region.
%
%Arguments
% s_point  Sequence of points on S^2, as a 2 x n_points array in spherical polar coordinates,
%          with longitude 0 <= s(1,p_idx) <= 2*pi, colatitude 0 <= s(2,p_idx) <= pi.
% N        Required number of regions, a positive integer.
%
%Examples
%
% >> points_s = eq_point_set_polar(2,8)
%
% points_s =
%
%          0    0.5236    1.5708    2.6180    3.6652    4.7124    5.7596         0
%          0    1.5708    1.5708    1.5708    1.5708    1.5708    1.5708    3.1416
%
% >> count_v = eq_count_points_by_s2_region(points_s, 8)
%
% count_v =
%
%      1     1     1     1     1     1     1     1
%
% >> count_v = eq_count_points_by_s2_region(points_s, 5)
%
% count_v =
%
%      1     2     2     2     1
%
% >> sum(count_v)
%
% ans =
%
%      8
%
% >> points_s = eq_point_set_polar(2,128);
%
% >> count_v = eq_count_points_by_s2_region(points_s, 8)
%
% count_v =
%
%     19    15    14    17    15    14    15    19
%
% >> sum(count_v)
%
% ans =
%
%    128
%
% >> count_v = eq_count_points_by_s2_region(points_s, 5)
%
% count_v =
%
%     19    29    32    29    19
%
% >> sum(count_v)
%
% ans =
%
%    128
%
%See also
% EQ_FIND_S2_REGION

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2012 Paul Leopardi
% $Revision 1.11 $ $Date 2012-01-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

r_idx = eq_find_s2_region(s_point, N);
count_v = histcounts(r_idx, 1:N+1);
