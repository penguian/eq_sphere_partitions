function r_idx = eq_find_s2_region(s_point, N)
%EQ_FIND_S2_REGION Given a set of points, partition S^2 into N regions
%and find the index of the region of containing each point
%
%Syntax
% r_idx = eq_find_s2_region(s_point, N)
%
%Description
% r_idx = eq_find_s2_region(s_point, N) does the following:
% 1) partitions the unit sphere S^2 into a sequence of N regions of
%    equal area and small diameter;
% 2) for each point in s_point, determines which region contains the point;
% 3) sets r_idx to be an array of length size(s_point, 2) containing the index of
%    the region corresponding to each point.
%
%Arguments
% s_point  Sequence of points on S^2, as a 2 x n_points array in spherical polar coordinates,
%          with longitude 0 <= s(1,p_idx) <= 2*pi, colatitude 0 <= s(2,p_idx) <= pi.
% N        Required number of regions, a positive integer.
%
%Examples
% > points_x=randn(3,8)
% points_x =
%
%     0.4889   -0.3034    0.8884   -0.8095    0.3252   -1.7115    0.3192   -0.0301
%     1.0347    0.2939   -1.1471   -2.9443   -0.7549   -0.1022    0.3129   -0.1649
%     0.7269   -0.7873   -1.0689    1.4384    1.3703   -0.2414   -0.8649    0.6277
%
% > points_s=cart2polar2(points_x)
% points_s =
%
%     1.1294    2.3722    5.3714    4.4441    5.1191    3.2013    0.7754    4.5321
%     1.0049    2.6491    2.2057    1.1306    0.5403    1.7107    2.6646    0.2609
%
% > r_idx = eq_find_s2_region(points_s,5)
% r_idx =
%
%      2     5     4     4     1     3     5     1
%
%See also
% EQ_COUNT_POINTS_BY_S2_REGION, LOOKUP_S2_REGION

% Copyright 2024 Paul Leopardi
% $Revision 1.12 $ $Date 2024-08-24 $
% Copyright 2012 Paul Leopardi
% $Revision 1.11 $ $Date 2012-01-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

[s_cap, n_regions] = eq_caps(2, N);
s_regions = eq_regions(2, N);
n_caps = length(s_cap);
c_regions = n_regions;
for c_idx = 2:n_caps
    c_regions(c_idx) = c_regions(c_idx - 1) + n_regions(c_idx);
end
r_idx = lookup_s2_region(s_point, s_regions, s_cap, c_regions);
