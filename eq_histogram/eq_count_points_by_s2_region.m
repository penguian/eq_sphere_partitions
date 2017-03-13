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
% N        Required number of regions, a positive integer.
% s_points Sequence of points on S^2, as a 2 x n_points array in spherical polar coordinates,
%          with longitude 0 <= s(1,p_idx) <= 2*pi, colatitude 0 <= s(2,p_idx) <= pi.
%
%Examples
%  > points_x=randn(3,8)
%  points_x =
%  
%    -1.500867  -0.994941  -0.466281  -2.419762  -0.027346  -0.362112   0.070776  -1.106204
%    -0.655976  -1.022458   0.669978  -1.066180  -0.025138   2.686040   1.195864  -0.087131
%    -0.273793  -0.965146  -1.040648   0.311548  -0.788428   1.790929   0.532415   0.919639
%  
%  > points_s=cart2polar2(points_x)
%  points_s =
%  
%     3.55364   3.94063   2.17881   3.55661   3.88495   1.70480   1.51168   3.22020
%     1.73642   2.16558   2.47645   1.45352   3.09452   0.98688   1.15258   0.87875
%  
%  > count_v=eq_count_points_by_s2_region(points_s,5)
%  count_v =
%  
%     1   2   3   0   2
%
%See also
% EQ_FIND_S2_REGION

% Copyright 2012 Paul Leopardi
% $Revision 1.11 $ $Date 2012-01-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

r_idx = eq_find_s2_region(s_point,N);
count_v = histc(r_idx,1:N);
