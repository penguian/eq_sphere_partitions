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
