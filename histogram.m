%Recursive Zonal Equal Area Sphere Partitioning: Histograms
%
% Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox.
% Release 1.12 2024-09-18
%
%Functions by category
%=====================
%
% Histogram utilities for S^2
%  eq_count_points_by_s2_region  Given a set of points, count points in each of N regions of S^2
%  eq_find_s2_region             Given a set of points, partition S^2 into N regions
%                                and find the index of the region of containing each point
%  in_s2_region                  Test that each of a set of points on S^2 is within a given region

% Copyright 2024 Paul Leopardi
% $Revision 1.12 $ $Date 2024-09-18 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.
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
%
% >> points_s = eq_point_set_polar(2,8)
%
% points_s =
%
%         0    0.5236    1.5708    2.6180    3.6652    4.7124    5.7596         0
%         0    1.5708    1.5708    1.5708    1.5708    1.5708    1.5708    3.1416
%
% >> r_idx = eq_find_s2_region(points_s,8)
%
% r_idx =
%
%      1     2     3     4     5     6     7     8
%
% >> r_idx = eq_find_s2_region(points_s,5)
%
% r_idx =
%
%      1     2     2     3     3     4     4     5
%
%See also
% EQ_COUNT_POINTS_BY_S2_REGION, LOOKUP_S2_REGION

% Copyright 2024 Paul Leopardi
% $Revision 1.12 $ $Date 2024-10-13 $
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
function result = in_s2_region(s_point,region)
%IN_S2_REGION Test that a point on S^2 is within a given region
%
%Syntax
% result = in_s2_region(s_point,region)
%
%Description
% result = in_s2_region(s_point,region)
% sets result to an array of length size(s_point,2) with each entry set to true (1)
% if the corresponding point in s_point is in region, false (0) otherwise.
%
%Arguments
% s_point  Sequence of points on S^2, as a 2 x n_points array in spherical polar coordinates,
%          with longitude 0 <= s(1,p_idx) <= 2*pi, colatitude 0 <= s(2,p_idx) <= pi.
% region   One region of S^2 as returned by eq_regions(2,N) for some positive integer N.
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
% >> s_regions = eq_regions(2,5);
% >> region = s_regions(:,:,3)
%
% region =
%     2.0944    4.1888
%     0.9273    2.2143
%
% >> result = in_s2_region(points_s, region)
%
% result =
%      0     0     0     1     1     0     0     0
%
%See also
%
% EQ_REGIONS, EQ_FIND_S2_REGION

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-09-17 $
% Copyright 2012 Paul Leopardi
% $Revision 1.11 $ $Date 2012-01-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_points = size(s_point, 2);
result = zeros(1, n_points);
for p_idx = 1:n_points
    result(p_idx) = false;
    longitude = s_point(1, p_idx);
    min_long = region(1, 1);
    max_long = region(1, 2);
    if min_long < longitude && longitude <= max_long
        result(p_idx) = true;
    else
        longitude = longitude + 2*pi;
        if min_long < longitude && longitude <= max_long
            result(p_idx) = true;
        end
    end
    colatitude = s_point(2, p_idx);
    min_colat = region(2, 1);
    max_colat = region(2, 2);
    if (min_colat == 0.0 && colatitude < min_colat) || ...
       (min_colat > 0.0 && colatitude <= min_colat) || ...
        max_colat < colatitude
        result(p_idx) = false;
    end
end
