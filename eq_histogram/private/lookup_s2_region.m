function r_idx = lookup_s2_region(s_point, s_regions, s_cap, c_regions)
%LOOKUP_S2_REGION For S^2, given sequences of points, regions, and cap colatitudes,
%find the index of the region containing each point.
%
%Syntax
% r_idx = lookup_s2_region(s_point,s_regions,s_cap,c_regions)
%
%Description
% r_idx = lookup_s2_region(s_point,s_regions,s_cap,c_regions) does the following:
% 1) for each point in s_point, determines which region of s_regions contains the point;
% 2) sets r_idx to be an array of length size(s_point,2) containing the index of
%    the region of s_regions corresponding to each point.
%
%Arguments
% s_point   Sequence of points on S^2, as a 2 x n_points array in spherical polar coordinates,
%           with longitude 0 <= s(1,p_idx) <= 2*pi, colatitude 0 <= s(2,p_idx) <= pi.
% s_regions Sequence of regions of S^2 as per eq_regions(2,N) where N == size(s_regions,2).
% s_cap     Sequence of cap colatitudes as per eq_caps(2,N) for the same N.
% c_regions Sequence of the cumulative number of regions of s_regions within each cap of c_cap.
%
%Examples
% > points_s = eq_point_set_polar(2,8)
% points_s =
%          0    0.5236    1.5708    2.6180    3.6652    4.7124    5.7596         0
%          0    1.5708    1.5708    1.5708    1.5708    1.5708    1.5708    3.1416
%
% > count_v = eq_count_points_by_s2_region(points_s, 8)
% count_v =
%      1     1     1     1     1     1     1     1
%
% > count_v = eq_count_points_by_s2_region(points_s, 5)
% count_v =
%      1     2     2     2     1
%
% > sum(count_v)
% ans =
%      8
%
% > points_s = eq_point_set_polar(2,128);
% > count_v = eq_count_points_by_s2_region(points_s, 8)
% count_v =
%     19    15    14    17    15    14    15    19
%
% > sum(count_v)
% ans =
%    128
%
% > count_v = eq_count_points_by_s2_region(points_s, 5)
% count_v =
%     19    29    32    29    19
%
% > sum(count_v)
% ans =
%    128
%
%
%See also
% EQ_REGIONS, EQ_CAPS, CUMSUM, LOOKUP_TABLE

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-09-07 $
% Copyright 2012 Paul Leopardi
% $Revision 1.11 $ $Date 2012-01-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_caps = length(s_cap);
if n_caps ~= length(c_regions)
    msg = 'LOOKUP_S2_REGION: Mismatch between length of s_cap (==%d) and length of c_regions (==%d)\n'
    fprintf(msg, n_caps, length(c_regions))
    r_idx = 0;
    return
end
n_regions = size(s_regions, 3);
if c_regions(n_caps) ~= n_regions
    msg = 'LOOKUP_S2_REGION: Mismatch between c_regions(end) (==%d) and length of s_regions (==%d)\n'
    fprintf(msg, c_regions(n_caps), n_regions)
    r_idx = 0;
    return
end
n_points = size(s_point, 2);
r_idx = zeros(1, n_points);
for p_idx = 1:n_points
    % Lookup by colatitude.
    c_idx = lookup_table(s_cap, s_point(2, p_idx));
    if c_idx > 0 && c_idx < n_caps - 1
        min_r_idx = c_regions(c_idx) + 1;
        max_r_idx = c_regions(c_idx + 1);
        s_longs = squeeze(s_regions(1, :, min_r_idx:max_r_idx));
        if s_longs(:, 1) >= 2*pi
            s_longs(:, 1) = s_longs(:, 1) - 2*pi;
        end
        n_longs = size(s_longs, 2);
	% Lookup by longitude.
        l_idx = mod(lookup_table(s_longs(2, :), s_point(1, p_idx)), n_longs);
        if s_point(1, p_idx) < s_longs(1, 1)
            l_idx = n_longs - 1;
        end
        r_idx(p_idx) = min_r_idx + l_idx;
    elseif c_idx == 0
        r_idx(p_idx) = 1;
    elseif c_idx >= n_caps - 1
        r_idx(p_idx) = n_regions;
    else
        r_idx(p_idx) = 0;
    end
end

