%Recursive Zonal Equal Area Sphere Partitioning: EQ partitions
%
% Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox.
% Release 1.12 2024-10-16
%
%Functions by category
%=====================
%
% Partitions
%  eq_caps                Partition a sphere into to nested spherical caps
%  eq_regions             Recursive zonal equal area (EQ) partition of sphere
%
% Point sets
%  eq_point_set           Center points of regions of EQ partition,
%                             in Cartesian coordinates
%  eq_point_set_polar     Center points of regions of an EQ partition,
%                             in spherical polar coordinates
%
% Partition options
%  partition_options      Options for EQ partition
%
% Illustration of algorithms
%  illustrate_eq_algorithm    Illustrate the EQ partition algorithm

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-16 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.
function [s_cap,n_regions] = eq_caps(dim,N)
%EQ_CAPS Partition a sphere into to nested spherical caps
%
%Syntax
% [s_cap,n_regions] = eq_caps(dim,N);
%
%Description
% [S_CAP,N_REGIONS] = EQ_CAPS(dim,N) does the following:
% 1) partitions the unit sphere S^dim into a list of spherical caps of
%    increasing colatitude and thus increasing area,
% 2) sets S_CAP to be an array of size (1 by N_COLLARS+2),
%    containing increasing colatitudes of caps, and
% 3) sets N_REGIONS to be an array of size (1 by N_COLLARS+2),
%    containing the integer number of regions in each corresponding zone of
%    S^dim.
%
% The argument N is assumed to be a positive integer.
%
%Notes
% The value N_COLLARS is a positive integer function of dim and N.
%
% S_CAP[1] is C_POLAR, the colatitude of the North polar cap.
% S_CAP[N_COLLARS+1] is pi-C_POLAR.
% S_CAP[N_COLLARS+2] is pi.
%
% N_REGIONS[1] is 1.
% N_REGIONS[N_COLLARS+2] is 1.
% The sum of N_REGIONS is N.
%
%Examples
%
% >> [s_cap,n_regions] = eq_caps(2,10)
%
% s_cap =
%
%     0.6435    1.5708    2.4981    3.1416
%
% n_regions =
%
%      1     4     4     1
%
% >> [s_cap,n_regions] = eq_caps(3,6)
%
% s_cap =
%
%     0.9845    2.1571    3.1416
%
% n_regions =
%
%      1     4     1
%
%See also
% EQ_REGIONS, EQ_POINT_SET_POLAR

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.1 $ $Date 2025-08-09 $
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
narginchk(2,2);
nargoutchk(2,2);
%
% dim is the number of dimensions
% N is the number of regions
%
if dim == 1
    %
    % We have a circle. Return the angles of N equal sectors.
    %
    sector = 1:N;
    %
    % Make dim==1 consistent with dim>1 by
    % returning the longitude of a sector enclosing the
    % cumulative sum of arc lengths given by summing n_regions.
    %
    s_cap = sector*2*pi/N;
    n_regions = ones(size(sector));
    %
elseif N == 1
    %
    % We have only one region, which must be the whole sphere.
    %
    s_cap = pi;
    n_regions = 1;
    %
else
    %
    % Given dim and N, determine c_polar,
    % the colatitude of the North polar spherical cap.
    %
    c_polar = polar_colat(dim,N);
    %
    % Given dim and N, determine the ideal angle for spherical collars.
    % Based on N, this ideal angle, and c_polar,
    % determine n_collars, the number of collars between the polar caps.
    %
    n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
    %
    % Given dim, N, c_polar and n_collars, determine r_regions,
    % a list of the ideal real number of regions in each collar,
    % plus the polar caps.
    % The number of elements is n_collars+2.
    % r_regions[1] is 1.
    % r_regions[n_collars+2] is 1.
    % The sum of r_regions is N.
    %
    r_regions = ideal_region_list(dim,N,c_polar,n_collars);
    %
    % Given N and r_regions, determine n_regions,
    % a list of the natural number of regions in each collar and
    % the polar caps.
    % This list is as close as possible to r_regions.
    % The number of elements is n_collars+2.
    % n_regions[1] is 1.
    % n_regions[n_collars+2] is 1.
    % The sum of n_regions is N.
    %
    n_regions = round_to_naturals(N,r_regions);
    %
    % Given dim, N, c_polar and n_regions, determine s_cap,
    % an increasing list of colatitudes of spherical caps which enclose the same area
    % as that given by the cumulative sum of regions.
    % The number of elements is n_collars+2.
    % s_cap[1] is c_polar.
    % s_cap[n_collars+1] is Pi-c_polar.
    % s_cap[n_collars+2] is Pi.
    %
    s_cap = cap_colats(dim,N,c_polar,n_regions);
    %
end
%
% end function
function points_x = eq_point_set(dim,N,varargin)
%EQ_POINT_SET Center points of regions of EQ partition, in Cartesian coordinates
%
%Syntax
% points_x = eq_point_set(dim,N,options);
%
%Description
% POINTS_X = EQ_POINT_SET(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
% partition S^dim (the unit sphere in dim+1 dimensional space) into N regions
% of equal area and small diameter, and
% 2) sets POINTS_X to be an array of size (dim+1 by N), containing the center
% points of each region.
% Each column of POINTS_X represents a point of S^dim, in Cartesian coordinates.
%
% The arguments dim and N must be positive integers.
%
% POINTS_X = EQ_POINT_SET(dim,N,'offset','extra') uses experimental extra offsets
% for S^2 and S^3 to try to minimize energy.
%
% POINTS_X = EQ_POINT_SET(dim,N,extra_offset) uses experimental extra offsets if
% extra_offset is true or non-zero.
%
%Notes
% Each region is defined as a product of intervals in spherical polar
% coordinates. The center point of a region is defined via the center points
% of each interval, with the exception of spherical caps and their descendants,
% where the center point is defined using the center of the spherical cap.
%
% If dim > 3, extra offsets are not used.
% For more details on options, see help partition_options.
%
%Examples
%
% >> points_x = eq_point_set(2,4)
%
% points_x =
%          0    0.0000   -0.0000    0.0000
%          0    1.0000   -1.0000         0
%     1.0000    0.0000    0.0000   -1.0000
%
% >> size(points_x)
%
% ans =
%
%      3     4
%
%See also
% PARTITION_OPTIONS, EQ_POINT_SET_POLAR, EQ_REGIONS, POLAR2CART

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-17 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
narginchk(2,4);
points_x = polar2cart(eq_point_set_polar(dim,N,varargin{:}));
%
% end function
function points_s = eq_point_set_polar(dim,N,varargin)
%EQ_POINT_SET_POLAR Center points of regions of an EQ partition
%
%Syntax
% points_s = eq_point_set_polar(dim,N,options);
%
%Description
% POINTS_S = EQ_POINT_SET_POLAR(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
% partition S^dim (the unit sphere in dim+1 dimensional space) into N regions
% of equal area and small diameter, and
% 2) sets POINTS_S to be an array of size (dim by N), containing the center
% points of each region. Each column of POINTS_S represents a point of S^dim,
% in spherical polar coordinates.
%
% The arguments dim and N must be positive integers.
%
% POINTS_S = EQ_POINT_SET_POLAR(dim,N,'offset','extra') uses experimental extra
% offsets for S^2 and S^3 to try to minimize energy. If dim > 3, extra offsets
% are not used.
%
% POINTS_S = EQ_POINT_SET_POLAR(dim,N,extra_offset) uses experimental extra
% offsets if extra_offset is true or non-zero.
%
%Notes
% Each region is defined as a product of intervals in spherical polar
% coordinates. The center point of a region is defined via the center points
% of each interval, with the exception of spherical caps and their descendants,
% where the center point is defined using the center of the spherical cap.
%
% For more details on options, see HELP PARTITION_OPTIONS.
%
%Examples
%
% >> points_s = eq_point_set_polar(2,4)
%
% points_s =
%
%          0    1.5708    4.7124         0
%          0    1.5708    1.5708    3.1416
%
% >> size(points_s)
%
% ans =
%
%      2     4
%
%See also
% PARTITION_OPTIONS, EQ_POINT_SET, EQ_REGIONS

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Function changed name from x2s2 to cart2polar2
% Optimize running time:
%   use slice assignments
%   trade space for time by using a cache
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
narginchk(2,4);
%
% dim is the number of dimensions
% N is the number of regions
%
% If the option 'offset' is 'extra', then use experimental extra offsets
% for S^2 and S^3 regions to try to minimize energy.
% The default is not to use extra offsets.
%
pdefault.extra_offset =  false;
if nargin < 3
    extra_offset = pdefault.extra_offset;
else
    popt = partition_options(pdefault, varargin{:});
    extra_offset = popt.extra_offset;
end
%
% Extra offset does not currently work for dim > 3,
% so quietly ignore this option in this case.
% Note that this also affects recursive calls to lower dimensions.
%
if dim > 3
    extra_offset = false;
end
%
% Check that dim and N are positive integers.
%
if ~( isnumeric(dim) && (dim >= 1) && (dim == floor(dim)) ) || ...
   ~( isnumeric(N) && (N >= 1) && (N == floor(N)) )
    fprintf('Usage: eq_point_set_polar(dim, N)\n');
    error('The arguments dim and N must be positive integers.');
end

if N == 1
    %
    % We have only one region, which must be the whole sphere.
    %
    points_s = zeros(dim,1);
    return;
end
%
% Start the partition of the sphere into N regions by partitioning
% to caps defined in the current dimension.
%
[a_cap, n_regions] = eq_caps(dim,N);
%
% a_cap is an increasing list of angles of the caps.
%
if dim == 1
    %
    % We have a circle and a_cap is an increasing list of angles of sectors,
    % with a_cap(k) being the cumulative arc length 2*pi/k.
    % The points are placed half way along each sector.
    %
    points_s = a_cap - pi/N;
    %
else
    %
    % We have a number of zones: two polar caps and a number of collars.
    % n_regions is the list of the number of regions in each zone.
    %
    n_collars = size(n_regions,2)-2;
    use_cache = dim >= 2;
    if use_cache
        cache_size = floor(n_collars/2);
        cache = cell(1,cache_size);
    end
    %
    % Start with the 'centre' point of the North polar cap.
    % This is the North pole.
    %
    points_s = zeros(dim,N);
    point_n = 2;
    %
    % Determine the 'centre' points for each collar.
    %
    if extra_offset && (dim == 3)
        R = eye(3);
    end
    if dim == 2
        offset = 0;
    end
    for collar_n = 1:n_collars
        %
        % a_top is the colatitude of the top of the current collar.
        %
        a_top = a_cap(collar_n);
        %
        % a_bot is the colatitude of the bottom of the current collar.
        %
        a_bot = a_cap(1+collar_n);
        %
        % n_in_collar is the number of regions in the current collar.
        %
        n_in_collar = n_regions(1+collar_n);
        %
        % The top and bottom of the collar are small (dim-1)-spheres,
        % which must be partitioned into n_in_collar regions.
        % Use eq_point_set_polar recursively to partition
        % the unit (dim-1)-sphere.
        % points_1 is the resulting list of points.
        %
        if use_cache
            twin_collar_n = n_collars-collar_n+1;
            if twin_collar_n <= cache_size && ...
                size(cache{twin_collar_n},2) == n_in_collar
                points_1 = cache{twin_collar_n};
            else
                points_1 = eq_point_set_polar(dim-1,n_in_collar,extra_offset);
                cache{collar_n} = points_1;
            end
        else
            points_1 = eq_point_set_polar(dim-1,n_in_collar,extra_offset);
        end
        %
        if extra_offset && (dim == 3) && (collar_n > 1)
            %
            % (Experimental)
            % Rotate 2-spheres to prevent alignment of north poles.
            %
            R = s2_offset(points_1)*R;
            points_1 = cart2polar2(R*polar2cart(points_1));
        end
        %
        % Given points_1, determine the 'centre' points for the collar.
        % Each point of points_1 is a 'centre' point on the (dim-1)-sphere.
        %
        % Angular 'centre' point;
        % The first angles are those of the current 'centre' point
        % of points_1, and the last angle in polar coordinates is the average of
        % the top and bottom angles of the collar,
        %
        a_point = (a_top+a_bot)/2;
        %
        point_1_n = 1:size(points_1,2);
        %
        if dim == 2
            %
            % The (dim-1)-sphere is a circle
            %
            points_s(1:dim-1,point_n+point_1_n-1) = mod(points_1(:,point_1_n)+2*pi*offset,2*pi);
            %
            % Given the number of sectors in the current collar and
            % in the next collar, calculate the next offset.
            % Accumulate the offset, and force it to be a number between 0 and 1.
            %
            offset = offset + circle_offset(n_in_collar,n_regions(2+collar_n),extra_offset);
            offset = offset - floor(offset);
        else
            points_s(1:dim-1,point_n+point_1_n-1) = points_1(:,point_1_n);
        end
        %
        points_s(dim, point_n+point_1_n-1) = a_point;
        point_n = point_n + size(points_1,2);
    end
    %
    % End with the 'centre' point of the bottom polar cap.
    %
    points_s(:,point_n) = zeros(dim,1);
    points_s(dim,point_n) = pi;
end
%
% end function
function [regions,dim_1_rot] = eq_regions(dim,N,varargin)
%EQ_REGIONS Recursive zonal equal area (EQ) partition of sphere
%
%Syntax
% [regions,dim_1_rot] = eq_regions(dim,N,options);
%
%Description
% REGIONS = EQ_REGIONS(dim,N) uses the recursive zonal equal area sphere
% partitioning algorithm to partition S^dim (the unit sphere in dim+1
% dimensional space) into N regions of equal area and small diameter.
%
% The arguments dim and N must be positive integers.
%
% The result REGIONS is a (dim by 2 by N) array, representing the regions
% of S^dim. Each element represents a pair of vertex points in spherical polar
% coordinates.
%
% Each region is defined as a product of intervals in spherical polar
% coordinates. The pair of vertex points regions(:,1,n) and regions(:,2,n) give
% the lower and upper limits of each interval.
%
% REGIONS = EQ_REGIONS(dim,N,'offset','extra') uses experimental extra
% offsets for S^2 and S^3 to try to minimize energy. If dim > 3, extra offsets
% are not used.
%
% REGIONS = EQ_REGIONS(dim,N,extra_offset) uses experimental extra offsets
% if extra_offset is true or non-zero.
%
% [REGIONS,DIM_1_ROT] = EQ_REGIONS(dim,N) also returns DIM_1_ROT, a cell
% array containing N rotation matrices, one per region, each of size dim by dim.
% These describe the R^dim rotation needed to place the region in its final
% position.
%
% [REGIONS,DIM_1_ROT] = EQ_REGIONS(dim,N,'offset','extra') partitions S^dim
% into N regions, using extra offsets, and also returning DIM_1_ROT, as above.
%
%Notes
% The output argument DIM_1_ROT is the only way to track the effect of the extra
% offset when dim == 3, because the R^3 rotation means that the boundary of a
% region generally no longer coincides with hyperplanes of colatitude and
% longitude. The function ILLUSTRATE_S3_PARTITION uses DIM_1_ROT.
%
% For more details on options, see HELP PARTITION_OPTIONS.
%
%Examples
%
% >> regions = eq_regions(2,4)
%
% regions(:,:,1) =
%          0    6.2832
%          0    1.0472
% regions(:,:,2) =
%          0    3.1416
%     1.0472    2.0944
% regions(:,:,3) =
%     3.1416    6.2832
%     1.0472    2.0944
% regions(:,:,4) =
%          0    6.2832
%     2.0944    3.1416
%
% >> size(regions)
%
% ans =
%      2     2     4
%
%See also
% PARTITION_OPTIONS, EQ_POINT_SET, EQ_POINT_SET_POLAR, PROJECT_S3_PARTITION

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Fix bug in assignment of dim_1_rot
% Documentation files renamed
% $Revision 1.02 $ $Date 2005-04-16 $
% Optimize running time:
%   move 'if nargout' blocks, refactor slice assignments
%   trade space for time by using a cache
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
narginchk(2,4);
nargoutchk(0,2);
%
% If the option 'offset' is 'extra', then use experimental extra offsets
% for S^2 and S^3 regions to try to minimize energy.
% The default is not to use extra offsets.
%
pdefault.extra_offset =  false;
if nargin < 3
    extra_offset = pdefault.extra_offset;
else
    popt = partition_options(pdefault, varargin{:});
    extra_offset = popt.extra_offset;
end
%
% Extra offset does not currently work for dim > 3,
% so quietly ignore this option in this case.
% Note that this also affects recursive calls to lower dimensions.
%
if dim > 3
    extra_offset = false;
end
%
% Check that dim and N are positive integers.
%
if ~( isnumeric(dim) && (dim >= 1) && (dim == floor(dim)) ) || ...
   ~( isnumeric(N) && (N >= 1) && (N == floor(N)) )
    fprintf('Usage: eq_regions(dim, N)\n');
    error('The arguments dim and N must be positive integers.');
end

if nargout > 1
    dim_1_rot = cell(1,N);
end

if N == 1
    %
    % We have only one region, which must be the whole sphere.
    %
    regions = zeros(dim,2,1);
    regions(:,:,1) = sphere_region(dim);
    if nargout > 1
        dim_1_rot{1} = eye(dim);
    end
    return;
end
%
% Start the partition of the sphere into N regions by partitioning
% to caps defined in the current dimension.
%
[s_cap, n_regions] = eq_caps(dim,N);
%
% s_cap is an increasing list of colatitudes of the caps.
%
if dim == 1
    %
    % We have a circle and s_cap is an increasing list of angles of sectors.
    %
    if nargout > 1
        R = eye(dim);
        for region_n = 1:N
            dim_1_rot{region_n} = R;
        end
    end
    %
    % Return a list of pairs of sector angles.
    %
    regions = zeros(dim,2,N);
    regions(:,1,2:N) = s_cap(1:N-1);
    regions(:,2,:)   = s_cap;
    %
else
    %
    % We have a number of zones: two polar caps and a number of collars.
    % n_regions is the list of the number of regions in each zone.
    %
    n_collars = size(n_regions,2)-2;
    use_cache = dim > 2;
    if use_cache
        cache_size = floor(n_collars/2);
        cache = cell(1,cache_size);
    end
    %
    % Start with the top cap.
    %
    regions = zeros(dim,2,N);
    regions(:,:,1) = top_cap_region(dim,s_cap(1));
    region_n = 1;
    %
    % Determine the dim-regions for each collar.
    %
    if (nargout > 1) || (extra_offset && (dim == 3))
        R = eye(dim);
    end
    if nargout > 1
        dim_1_rot{1} = R;
    end
    if dim == 2
        offset = 0;
    end
    for collar_n = 1:n_collars
        %
        % c_top is the colatitude of the top of the current collar.
        %
        c_top = s_cap(collar_n);
        %
        % c_bot is the colatitude of the bottom of the current collar.
        %
        c_bot = s_cap(1+collar_n);
        %
        % n_in_collar is the number of regions in the current collar.
        %
        n_in_collar = n_regions(1+collar_n);
        %
        % The top and bottom of the collar are small (dim-1)-spheres,
        % which must be partitioned into n_in_collar regions.
        % Use eq_regions recursively to partition
        % the unit (dim-1)-sphere.
        % regions_1 is the resulting list of (dim-1)-region pairs.
        %
        if use_cache
            twin_collar_n = n_collars-collar_n+1;
            if twin_collar_n <= cache_size && ...
                size(cache{twin_collar_n},3) == n_in_collar
                regions_1 = cache{twin_collar_n};
            else
                regions_1 = eq_regions(dim-1,n_in_collar,extra_offset);
                cache{collar_n} = regions_1;
            end
        else
            regions_1 = eq_regions(dim-1,n_in_collar,extra_offset);
        end
        %
        if extra_offset && (dim == 3) && (collar_n > 1)
            %
            % (Experimental)
            % Rotate 2-spheres to prevent alignment of north poles.
            %
            R = s2_offset(centres_of_regions(regions_1))*R;
        end
        %
        % Given regions_1, determine the dim-regions for the collar.
        % Each element of regions_1 is a (dim-1)-region pair for
        % the (dim-1)-sphere.
        %
        if nargout > 1
            for region_1_n = 1:size(regions_1,3)
                dim_1_rot{region_n+region_1_n} = R;
            end
        end
        if dim == 2
            %
            % The (dim-1)-sphere is a circle
            % Offset each sector angle by an amount which accumulates over
            % each collar.
            %
            for region_1_n = 1:size(regions_1,3)
                %
                % Top of 2-region;
                % The first angle is the longitude of the top of
                % the current sector of regions_1, and
                % the second angle is the top colatitude of the collar.
                %
                r_top = [mod(regions_1(1,1,region_1_n)+2*pi*offset,2*pi); c_top];
                %
                % Bottom of 2-region;
                % The first angle is the longitude of the bottom of
                % the current sector of regions_1, and
                % the second angle is the bottom colatitude of the collar.
                %
                r_bot = [mod(regions_1(1,2,region_1_n)+2*pi*offset,2*pi); c_bot];
                if r_bot(1) < r_top(1)
                   r_bot(1) = r_bot(1) + 2*pi;
                end
                region_n = region_n+1;
                regions(:,:,region_n) = [r_top,r_bot];
            end
            %
            % Given the number of sectors in the current collar and
            % in the next collar, calculate the next offset.
            % Accumulate the offset, and force it to be a number between 0 and 1.
            %
            offset = offset + circle_offset(n_in_collar,n_regions(2+collar_n),extra_offset);
            offset = offset - floor(offset);
        else
            for region_1_n = 1:size(regions_1,3)
                %
                region_n = region_n+1;
                %
                % Dim-region;
                % The first angles are those of the current (dim-1) region of regions_1.
                %
                regions(1:dim-1,:,region_n) = regions_1(:,:,region_1_n);
                %
                % The last angles are the top and bottom colatitudes of the collar.
                %
                regions(dim,:,region_n) = [c_top,c_bot];
            end
        end
    end
    %
    % End with the bottom cap.
    %
    regions(:,:,N) = bot_cap_region(dim,s_cap(1));
    if nargout > 1
        dim_1_rot{N} = eye(dim);
    end
end
%
% end function
function illustrate_eq_algorithm(dim,N,varargin)
%ILLUSTRATE_EQ_ALGORITHM Illustrate the EQ partition algorithm
%
%Syntax
% illustrate_eq_algorithm(dim,N,options);
%
%Description
% ILLUSTRATE_EQ_ALGORITHM(dim,N) illustrates the recursive zonal equal area
% sphere partitioning algorithm, which partitions S^dim (the unit sphere in
% dim+1 dimensional space) into N regions of equal area and small diameter.
%
% The illustration consists of four subplots:
% 1. Steps 1 and 2
% 2. Steps 3 to 5
% 3. Steps 6 and 7
% 4. Lower dimensional partitions (if dim == 2 or dim == 3)
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'offset','extra') uses experimental extra
% offsets for S^2 and S^3. If dim > 3, extra offsets are not used.
% For more detail on partition options, see HELP PARTITION_OPTIONS.
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,options) also recognizes a number of
% illustration options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'fontsize',size)
% Font size used in titles (numeric, default 16).
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'title','long')
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'title','short')
% Use long or short titles (default 'short').
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'proj','stereo')
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'proj','eqarea')
% Use stereographic or equal area projection (default 'stereo').
%
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'points','show')
% ILLUSTRATE_EQ_ALGORITHM(dim,N,'points','hide')
% Show or hide center points (default 'show').
%
% See examples below.
% For more detail on illustration options, see HELP ILLUSTRATION_OPTIONS.
%
%Notes
% The step numbers refer to the following steps of the the recursive zonal
% equal area sphere partitioning algorithm, which partition the sphere into
% zones.
%
% 1. Determine the colatitudes of the North and South polar caps.
% 2. Determine an ideal collar angle.
% 3. Use the angle between the North and South polar caps and the ideal collar
%    angle to determine an ideal number of collars.
% 4. Use a rounding procedure to determine the actual number of collars,
%    given the ideal number of collars.
% 5. Create a list containing the ideal number of regions in each collar.
% 6. Use a rounding procedure to create a list containing the actual number of
%    regions in each collar, given the list containing the ideal number of
%    regions.
% 7. Create a list containing the colatitude of the top of each zone,
%    given the list containing the actual number of regions in each collar,
%    and the colatitudes of the polar caps.
%
%Examples
%
% >> illustrate_eq_algorithm(3,99)
% >> illustrate_eq_algorithm(3,99,'offset','extra','proj','eqarea')
% >> illustrate_eq_algorithm(3,99,'proj','eqarea','points','hide')
%
%See also
% PARTITION_OPTIONS, ILLUSTRATION_OPTIONS, SUBPLOT

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.1 $ $Date 2025-08-09 $
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

pdefault.extra_offset =  false;

popt = partition_options(pdefault, varargin{:});

gdefault.fontsize = 16;
gdefault.show_title =    true;
gdefault.long_title =    false;
gdefault.stereo =        false;
gdefault.show_points =   true;

gopt = illustration_options(gdefault, varargin{:});
opt_args = option_arguments(popt,gopt);

subplot(2,2,1);axis off
illustrate_steps_1_2(dim,N,opt_args);

subplot(2,2,2);axis off
illustrate_steps_3_5(dim,N,opt_args);

subplot(2,2,3);axis off
illustrate_steps_6_7(dim,N,opt_args);

subplot(2,2,4);axis off
cla

gopt.fontsize = 32;
switch dim
case 2
    opt_args = option_arguments(popt,gopt);
    project_s2_partition(N,opt_args{:});
case 3
    opt_args = option_arguments(popt,gopt);
    [~,m] = eq_caps(dim,N);
    max_collar = min(4,size(m,2)-2);
    for k = 1:max_collar
        subn = 9+2*k-mod(k-1,2);
        subplot(4,4,subn);axis off
        project_s2_partition(m(1+k),opt_args{:});
    end
end
%
% end function

function illustrate_steps_1_2(dim,N,varargin)
% Illustrate steps 1 and 2 of the EQ partition of S^dim into N regions;
%
% illustrate_steps_1_2(dim,N,options);

gdefault.fontsize = 14;
gdefault.show_title =    true;
gdefault.long_title =    false;

gopt = illustration_options(gdefault, varargin{:});
h = 0:1/90:1;
% Plot a circle to represent dth coordinate of S^d
Phi = h*2*pi;
plot(sin(Phi),cos(Phi),'k','LineWidth',1)
axis equal;axis off;hold on

c_polar = polar_colat(dim,N);

k = -1:1/20:1;
j = ones(size(k));

% Plot the bounding parallels of the polar caps
plot(sin(c_polar)*k, cos(c_polar)*j,'r','LineWidth',2)
plot(sin(c_polar)*k,-cos(c_polar)*j,'r','LineWidth',2)

% Plot the North-South axis
plot(zeros(size(j)),k,'b','LineWidth',1)
% Plot the polar angle
plot(sin(c_polar)*h,cos(c_polar)*h,'b','LineWidth',2)

text(0.05,2/3,'\theta_c','Fontsize',gopt.fontsize);

% Plot the ideal collar angle
Delta_I = ideal_collar_angle(dim,N);
theta = c_polar + Delta_I;
plot(sin(theta)*h,cos(theta)*h,'b','LineWidth',2)

mid = c_polar + Delta_I/2;
text(sin(mid)*2/3,cos(mid)*2/3,'\Delta_I','Fontsize',gopt.fontsize);

% Plot an arc to indicate angles
theta = h*(c_polar + Delta_I);
plot(sin(theta)/5,cos(theta)/5,'b','LineWidth',1)

text(-0.9,-0.1,sprintf('V(\\theta_c) = V_R \n    = \\sigma(S^{%d})/%d',dim,N),...
    'Fontsize',gopt.fontsize);

caption_angle = min(mid + 2*Delta_I,pi-c_polar);
text(sin(caption_angle)/3,cos(caption_angle)/3,sprintf('\\Delta_I = V_R^{1/%d}',dim),...
    'Fontsize',gopt.fontsize);

if gopt.show_title
    title_str = sprintf('EQ(%d,%d) Steps 1 to 2\n',dim,N);
    title(title_str,'Fontsize',gopt.fontsize);
end

hold off
%
% end function

function illustrate_steps_3_5(dim,N,varargin)
% Illustrate steps 3 to 5 of the EQ partition of S^dim into N regions;
%
% illustrate_steps_3_5(dim,N,options);

gdefault.fontsize = 14;
gdefault.show_title =    true;
gdefault.long_title =    false;

gopt = illustration_options(gdefault, varargin{:});

h = 0:1/90:1;
Phi = h*2*pi;
plot(sin(Phi),cos(Phi),'k','LineWidth',1)
axis equal;axis off;hold on

c_polar = polar_colat(dim,N);
n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
r_regions = ideal_region_list(dim,N,c_polar,n_collars);
s_cap = cap_colats(dim,N,c_polar,r_regions);

k = -1:1/20:1;
j = ones(size(k));
plot(sin(c_polar)*k, cos(c_polar)*j,'r','LineWidth',2);

plot(zeros(size(j)),k,'b','LineWidth',1)

for collar_n = 0:n_collars
    zone_n = 1+collar_n;
    theta = s_cap(zone_n);
    plot(sin(theta)*h,cos(theta)*h,'b','LineWidth',2);
    theta_str = sprintf('\\theta_{F,%d}',zone_n);
    text(sin(theta)*1.1,cos(theta)*1.1,theta_str,'Fontsize',gopt.fontsize);
    if collar_n ~= 0
        plot(sin(theta)*k, cos(theta)*j,'r','LineWidth',2);
        theta_p = s_cap(collar_n);
        arc = theta_p + (theta-theta_p)*h;
        plot(sin(arc)/5,cos(arc)/5,'b','LineWidth',1);
        mid = (theta_p + theta)/2;
        text(sin(mid)/2,cos(mid)/2,'\Delta_F','Fontsize',gopt.fontsize);
        y_str = sprintf('y_{%d} = %3.1f...',collar_n,r_regions(zone_n));
        text(-sin(mid)+1/20,cos(mid)+(mid-pi)/30,y_str,'Fontsize',gopt.fontsize);
    end
end
if gopt.show_title
    title_str = sprintf('EQ(%d,%d) Steps 3 to 5\n',dim,N);
    title(title_str,'Fontsize',gopt.fontsize);
end
hold off
%
% end function

function illustrate_steps_6_7(dim,N,varargin)
% Illustrate steps 6 to 7 of the EQ partition of S^dim into N regions;
%
% illustrate_steps_6_7(dim,N,options);

gdefault.fontsize = 14;
gdefault.show_title =    true;
gdefault.long_title =    false;

gopt = illustration_options(gdefault, varargin{:});

h = 0:1/90:1;
Phi = h*2*pi;
plot(sin(Phi),cos(Phi),'k','LineWidth',1)
axis equal;axis off;hold on

c_polar = polar_colat(dim,N);
n_collars = num_collars(N,c_polar,ideal_collar_angle(dim,N));
r_regions = ideal_region_list(dim,N,c_polar,n_collars);
n_regions = round_to_naturals(N,r_regions);
s_cap = cap_colats(dim,N,c_polar,n_regions);

k = -1:1/20:1;
j = ones(size(k));
plot(sin(c_polar)*k, cos(c_polar)*j,'r','LineWidth',2);

plot(zeros(size(j)),k,'b','LineWidth',1)

for collar_n = 0:n_collars
    zone_n = 1+collar_n;
    theta = s_cap(zone_n);
    plot(sin(theta)*h,cos(theta)*h,'b','LineWidth',2);
    theta_str = sprintf('\\theta_{%d}',zone_n);
    text(sin(theta)*1.1,cos(theta)*1.1,theta_str,'Fontsize',gopt.fontsize);
    if collar_n ~= 0
        plot(sin(theta)*k, cos(theta)*j,'r','LineWidth',2);
        theta_p = s_cap(collar_n);
        arc = theta_p + (theta-theta_p)*h;
        plot(sin(arc)/5,cos(arc)/5,'b','LineWidth',1);
        mid = (theta_p + theta)/2;
        Delta_str = sprintf('\\Delta_{%i}',collar_n);
        text(sin(mid)/2,cos(mid)/2,Delta_str,'Fontsize',gopt.fontsize);
        m_str = sprintf('m_{%d} =%3.0f',collar_n,n_regions(zone_n));
        text(-sin(mid)+1/20,cos(mid)+(mid-pi)/30,m_str,'Fontsize',gopt.fontsize);
    end
end
if gopt.show_title
    title_str = sprintf('EQ(%d,%d) Steps 6 to 7\n',dim,N);
    title(title_str,'Fontsize',gopt.fontsize);
end
hold off
%
% end function

function arg = option_arguments(popt,gopt)

k = 1;
if isfield(popt,'extra_offset')
    arg{k} = 'offset';
    if popt.extra_offset
        arg{k+1} = 'extra';
    else
        arg{k+1} = 'normal';
    end
    k = k+2;
end

if isfield(gopt,'fontsize')
    arg{k} = 'fontsize';
    arg{k+1} = gopt.fontsize;
    k = k+2;
end

if isfield(gopt,'stereo')
    arg{k} = 'proj';
    if gopt.stereo
        arg{k+1} = 'stereo';
    else
        arg{k+1} = 'eqarea';
    end
    k = k+2;
end

if isfield(gopt,'show_title')
    arg{k} = 'title';
    if gopt.show_title
        if isfield(gopt,'long_title')
            if gopt.long_title
                arg{k+1} = 'long';
            else
                arg{k+1} = 'short';
            end
        else
            arg{k+1} = 'show';
        end
    else
        arg{k+1} = 'none';
    end
    k = k+2;
elseif isfield(gopt,'long_title')
    arg{k} = 'title';
    if gopt.long_title
        arg{k+1} = 'long';
    else
        arg{k+1} = 'short';
    end
    k = k+2;
end


if isfield(gopt,'show_points')
    arg{k} = 'points';
    if gopt.show_points
        arg{k+1} = 'show';
    else
        arg{k+1} = 'hide';
    end
    k = k+2;
end

if isfield(gopt,'show_surfaces')
    arg{k} = 'surf';
    if gopt.show_surfaces
        arg{k+1} = 'show';
    else
        arg{k+1} = 'hide';
    end
end

function popt = partition_options(pdefault, varargin)
%PARTITION_OPTIONS Options for EQ partition
%
%Syntax
% popt = partition_options(pdefault,options);
%
%Description
% POPT = PARTITION_OPTIONS(PDEFAULT,options) collects partition options,
% specified as name, value pairs, and places these into the structure POPT.
% The structure PDEFAULT is used to define default option values.
%
% The structures pdefault and popt may contain the following fields:
% extra_offset:  boolean
%
% The following partition options are available.
%
% 'offset':      Control extra rotation offsets for S^2 and S^3 regions.
%     'extra':   Use extra rotation offsets for S^2 and S^3 regions, to try
%                to minimize energy.
%                Sets opt.extra_offset to true.
%     'normal':  Do not use extra offsets
%                Sets opt.extra_offset to false.
%
% Some shortcuts are also provided.
% POPT = PARTITION_OPTIONS(pdefault) just sets POPT to PDEFAULT.
%
% The following are equivalent to PARTITION_OPTIONS(PDEFAULT,'offset','extra'):
% PARTITION_OPTIONS(PDEFAULT,true)
% PARTITION_OPTIONS(PDEFAULT,'extra')
%
% The following are equivalent to PARTITION_OPTIONS(PDEFAULT,'offset','normal'):
% PARTITION_OPTIONS(PDEFAULT,false)
% PARTITION_OPTIONS(PDEFAULT,'normal')
%
%Examples
%
% >> pdefault.extra_offset = false;
% >> popt = partition_options(pdefault,'offset','extra')
%
% popt =
%
%     extra_offset: 1
%
% >> popt = partition_options(pdefault,false)
%
% popt =
%     extra_offset: 0

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

popt = pdefault;
nargs = length(varargin);

if nargs == 1
    %
    % Short circuit: single argument is value of extra_offset
    %
    value = varargin{1};
    switch value
    case true
        popt.extra_offset = true;
    case false
        popt.extra_offset = false;
    case 'extra'
        popt.extra_offset = true;
    case 'normal'
        popt.extra_offset = false;
    otherwise
        value_error(value,varargin{:});
    end
    return;
end

nopts = floor(nargs/2);
opt_args = varargin(1:2:2*nopts-1);
for k=1:nopts
    if ~ischar([opt_args{k}])
        fprintf('Option names must be character strings\n');
        option_error(varargin{:});
    end
end
opt_vals = varargin(2:2:2*nopts);

option_name = 'offset';
pos = find(strcmp(option_name,opt_args));
if ~isempty(pos)
    if (isscalar(pos))
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'extra'
        popt.extra_offset = true;
    case 'normal'
        popt.extra_offset = false;
    otherwise
        value_error(value,varargin{:});
    end
end

function duplicate_error(option_name,varargin)
fprintf('Duplicate option %s\n',option_name);
option_error(varargin{:});
%
% end function

function value_error(value,varargin)
fprintf('Invalid option value ');
disp(value);
option_error(varargin{:});
%
% end function

function option_error(varargin)
fprintf('Error in options:\n');
disp(varargin);
error('Please check "help partition_options" for options');
%
% end function
