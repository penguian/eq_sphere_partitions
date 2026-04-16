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
% s_regions Sequence of regions of S^2 as per eq_regions(2,N) where N == size(s_regions,3).
% s_cap     Sequence of cap colatitudes as per eq_caps(2,N) for the same N.
% c_regions Sequence of the cumulative number of regions of s_regions within each cap of s_cap.
%
%Examples
% > cd eq_histogram/private
% > points_s = eq_point_set_polar(2,8)
% points_s =
%          0    0.5236    1.5708    2.6180    3.6652    4.7124    5.7596         0
%          0    1.5708    1.5708    1.5708    1.5708    1.5708    1.5708    3.1416
%
% > N = 8;
% > s_regions = eq_regions(2, N);
% > [s_cap, n_regions] = eq_caps(2, N);
% > c_regions = cumsum(n_regions);
% > r_idx = lookup_s2_region(points_s, s_regions, s_cap, c_regions)
% r_idx =
%      1     2     3     4     5     6     7     8
%
% > N = 5;
% > s_regions = eq_regions(2, N);
% > [s_cap, n_regions] = eq_caps(2, N);
% > c_regions = cumsum(n_regions);
% > r_idx = lookup_s2_region(points_s, s_regions, s_cap, c_regions)
% r_idx =
%      1     2     2     3     3     4     4     5
%
%See also
% EQ_REGIONS, EQ_CAPS, CUMSUM, LOOKUP_TABLE

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.3 $ $Date 2026-04-14 $
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-09-18 $
% Copyright 2012 Paul Leopardi
% $Revision 1.11 $ $Date 2012-01-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_caps = length(s_cap);
if n_caps ~= length(c_regions)
    msg = 'LOOKUP_S2_REGION: Mismatch between length of s_cap (==%d) and length of c_regions (==%d)\n';
    error(msg, n_caps, length(c_regions))
end
n_regions = size(s_regions, 3);
if c_regions(n_caps) ~= n_regions
    msg = 'LOOKUP_S2_REGION: Mismatch between c_regions(end) (==%d) and length of s_regions (==%d)\n';
    error(msg, c_regions(n_caps), n_regions)
end

tol = eps * 2^5;
n_points = size(s_point, 2);
r_idx = zeros(1, n_points);

if n_points == 0
    return
end

for p_idx = 1:n_points
    lat = s_point(2, p_idx);
    % Find cap index: first cap boundary >= lat - tol
    c_idx = find(s_cap >= lat - tol, 1);
    if isempty(c_idx)
        c_idx = n_caps;
    end

    % min_r_idx: index of first region in this cap band
    if c_idx == 1
        min_r_idx = 1;
        n_longs = c_regions(1);
    else
        min_r_idx = c_regions(c_idx-1) + 1;
        n_longs = c_regions(c_idx) - (min_r_idx - 1);
    end

    if n_longs > 1
        start_off = min_r_idx - 1;
        % s_regions is (2, 2, N). Longitude is dimension 1.
        % regions(1, 1, :) are starts, regions(1, 2, :) are ends.
        phi0 = s_regions(1, 1, start_off + 1);
        ends = squeeze(s_regions(1, 2, start_off + 1 : start_off + n_longs))';

        % Translate point longitude to [0, 2*pi) relative to first region start
        pts_long_translated = mod(s_point(1, p_idx) - phi0, 2*pi);
        % Boundary policy: longitude 0 (relative) is treated as 2*pi
        % to fall into the last region of the collar.
        if pts_long_translated <= tol
            pts_long_translated = 2*pi;
        end

        % Translate boundaries (ends) to [0, 2*pi)
        ends_translated = mod(ends - phi0, 2*pi);

        % The ends are monotonically increasing in [0, 2*pi) EXCEPT for the last one
        % which wraps exactly to 0. Set it to 2*pi for monotonicity.
        if ends_translated(end) <= tol
            ends_translated(end) = 2*pi;
        end

        % Find sector index in collar: first end >= pts_long_translated - tol
        l_idx = find(ends_translated >= pts_long_translated - tol, 1);

        % Overshoot protection
        if isempty(l_idx) || l_idx > n_longs
            l_idx = 1;
        end

        r_idx(p_idx) = min_r_idx + l_idx - 1;
    else
        % Pole cap logic
        r_idx(p_idx) = min_r_idx;
    end
end
