%Recursive Zonal Equal Area Sphere Partitioning: Properties of EQ partitions
%
% Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox.
% Release 1.12 2024-10-13
%
%Functions by category
%=====================
%
% Area
%  eq_area_error          Total area error and max area error per region of an
%                             EQ partition
%
% Diameter
%  eq_diam_bound          Maximum per-region diameter bound of EQ partition
%  eq_vertex_diam         Maximum vertex diameter of EQ partition
%
%  eq_diam_coeff          Coefficients of diameter bound and vertex diameter of
%                             EQ partition
%  eq_vertex_diam_coeff   Coefficient of maximum vertex diameter of EQ partition
%
% Hook for user-defined properties
%  eq_regions_property    Property of regions of an EQ partition

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
function [total_error, max_error] = eq_area_error(dim,N)
%EQ_AREA_ERROR Total area error and max area error per region of an EQ partition
%
%Syntax
% [total_error, max_error] = eq_area_error(dim,N)
%
%Description
% [TOTAL_ERROR, MAX_ERROR] = EQ_AREA_ERROR(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) sets TOTAL_ERROR to be the absolute difference between the total area of
%    all regions of the partition, and the area of S^dim, and
% 3) sets MAX_ERROR to be the maximum absolute difference between the area of
%    any region of the partition, and the ideal area of a region as given by
%    AREA_OF_IDEAL_REGION(dim,N), which is 1/N times the area of S^dim.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The results TOTAL_ERROR and MAX_ERROR will be arrays of the same size as N.
%
%Examples
%
% >> [total_error, max_error] = eq_area_error(2,10)
%
% total_error =
%
%    1.7764e-15
%
% max_error =
%
%    4.4409e-16
%
% >> [total_error, max_error] = eq_area_error(3,1:6)
%
% total_error =
%
%    1.0e-12 *
%     0.0036    0.0036    0.1847    0.0142    0.0142    0.2132
%
% max_error =
%
%    1.0e-12 *
%     0.0036    0.0018    0.1954    0.0284    0.0440    0.0777
%
%See also
% EQ_REGIONS, AREA_OF_SPHERE, AREA_OF_IDEAL_REGION

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
% Flatten N into a row vector.
%
shape = size(N);
n_partitions = prod(shape);
N = reshape(N,1,n_partitions);

total_error = zeros(size(N));
max_error   = zeros(size(N));
sphere_area = area_of_sphere(dim);

for partition_n = 1:n_partitions
    n = N(partition_n);
    regions = eq_regions(dim,n);
    ideal_area = area_of_ideal_region(dim,n);
    total_area = 0;
    for region_n = 1:size(regions,3)
        area = area_of_region(regions(:,:,region_n));
        total_area = total_area + area;
        region_error = abs(area - ideal_area);
        if region_error > max_error(partition_n)
            max_error(partition_n) = region_error;
        end
    end
    total_error(partition_n) = abs(sphere_area - total_area);
end
%
% Reshape output to same array size as original N.
%
total_error = reshape(total_error,shape);
max_error = reshape(max_error,shape);
%
% end function

function area = area_of_region(region)
%AREA_OF_REGION Area of given region
%
% area = area_of_region(region);

dim = size(region,1);
s_top = region(dim,1);
s_bot = region(dim,2);
if dim > 1
    area = area_of_collar(dim, s_top, s_bot)*area_of_region(region(1:dim-1,:))/area_of_sphere(dim-1);
else
    if s_bot == 0
        s_bot = 2*pi;
    end
    if s_top == s_bot
        s_bot = s_top + 2*pi;
    end
    area = s_bot - s_top;
end
%
% end function
function diam_bound = eq_diam_bound(dim,N)
%EQ_DIAM_BOUND Maximum per-region diameter bound of EQ partition
%
%Syntax
% diam_bound = eq_diam_bound(dim,N);
%
%Description
% DIAM_BOUND = EQ_DIAM_BOUND(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) sets DIAM_BOUND to be the maximum of the per-region diameter bound over
%    all the regions of the partition.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result DIAM_BOUND will an array of the same size as N.
%
%Examples
%
% >> diam_bound = eq_diam_bound(2,10)
%
%  diam_bound =
%
%      1.6733
%
% >> diam_bound = eq_diam_bound(3,1:6)
%
%  diam_bound =
%
%       2     2     2     2     2     2
%
%See also
% EQ_VERTEX_DIAM, EQ_DIAM_COEFF, EQ_REGIONS_PROPERTY

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

diam_bound = eq_regions_property(@max_diam_bound_of_regions,dim,N);
%
% end function
function [bound_coeff,vertex_coeff] =  eq_diam_coeff(dim,N)
%EQ_DIAM_COEFF Coefficients of diameter bound and vertex diameter of EQ partition
%
%Syntax
% [bound_coeff,vertex_coeff] = eq_diam_coeff(dim,N);
%
%Description
% [BOUND_COEFF,VERTEX_COEFF] = EQ_DIAM_COEFF(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) finds the maximum of the per-region diameter bound over all the regions
%    of the partition,
% 3) sets BOUND_COEFF to be the diameter bound coefficient, defined as the
%    solution to
%
%    max_diam_bound == BOUND_COEFF N^(-1/dim),
%
% 4) optionally finds the maximum vertex diameter over all the regions of the
%    partition, and
% 5) optionally sets VERTEX_COEFF to be the vertex diameter coefficient,
%    defined as the solution to
%
%    max_vertex_diam == VERTEX_COEFF N^(-1/dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result BOUND_COEFF and the optional result VERTEX_COEFF will be arrays of
% the same size as N.
%
%Examples
%
% >> bound_coeff = eq_diam_coeff(2,10)
%
%  bound_coeff =
%
%      5.2915
%
% >> [bound_coeff,vertex_coeff]=eq_diam_coeff(3,1:6)
%
%  bound_coeff =
%
%      2.0000    2.5198    2.8845    3.1748    3.4200    3.6342
%
%  vertex_coeff =
%
%      2.0000    2.5198    2.8845    3.1748    3.4200    3.6342
%
%See also
% EQ_DIAM_BOUND, EQ_VERTEX_DIAM, EQ_REGIONS, EQ_VERTEX_DIAM_COEFF

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of NSW.
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
nargoutchk(0,2);

if nargout < 2
    bound_coeff = eq_diam_bound(dim,N) .* N.^(1/dim);
else
    %
    % Flatten N into a row vector.
    %
    shape = size(N);
    n_partitions = prod(shape);
    N = reshape(N,1,n_partitions);

    bound_coeff =  zeros(size(N));
    vertex_coeff = zeros(size(N));
    for partition_n = 1:n_partitions
        n = N(partition_n);
        regions = eq_regions(dim,n);
        scale = n^(1/dim);
        bound_coeff(partition_n) =  max_diam_bound_of_regions(regions)  * scale;
        vertex_coeff(partition_n) = max_vertex_diam_of_regions(regions) * scale;
    end
    %
    % Reshape output to same array size as original N.
    %
    bound_coeff =  reshape(bound_coeff,shape);
    vertex_coeff = reshape(vertex_coeff,shape);
end
%
%end function
function property = eq_regions_property(fhandle,dim,N)
%EQ_REGIONS_PROPERTY Property of regions of an EQ partition
%
%Syntax
% property = eq_regions_property(fhandle,dim,N);
%
%Description
% PROPERTY = EQ_REGIONS_PROPERTY(FHANDLE,dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) calls the function defined by FHANDLE which is expected to use the
%    regions to calculate the value of the result, PROPERTY.
%
% The argument FHANDLE must be a function handle. The function specified by
% FHANDLE must take as its argument a single (dim by 2 by N) array,
% representing N regions of S^dim, in spherical polar coordinates, and must
% return a single value based on these regions.
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result PROPERTY will be an array of the same size as N.
%
%Examples
%
% See code in Matlab M files eq_diam_bound.m, eq_vertex_diam.m.
%
%See also
% EQ_REGIONS, FEVAL, EQ_DIAM_BOUND, EQ_VERTEX_DIAM

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
% Flatten N into a row vector.
%
shape = size(N);
n_partitions = prod(shape);
N = reshape(N,1,n_partitions);

property = zeros(size(N));
for partition_n = 1:n_partitions
    regions = eq_regions(dim,N(partition_n));
    property(partition_n) = feval(fhandle,regions);
end
%
% Reshape output to same array size as original N.
%
property = reshape(property,shape);
%
% end function
function coeff =  eq_vertex_diam_coeff(dim,N)
%EQ_VERTEX_DIAM_COEFF Coefficient of maximum vertex diameter of EQ partition
%
%Syntax
% coeff = eq_vertex_diam_coeff(dim,N);
%
%Description
% COEFF = EQ_VERTEX_DIAM_COEFF(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) finds the maximum vertex diameter over all the regions of the partition,
% 3) sets COEFF to be the vertex diameter coefficient,
%    defined as the solution to
%
%    max_vertex_diam == COEFF N^(-1/dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result COEFF will an array of the same size as N.
%
%Examples
%
% >> coeff = eq_vertex_diam_coeff(2,10)
%
%  coeff =
%
%      4.4721
%
% >> coeff = eq_vertex_diam_coeff(3,1:6)
%
%  coeff =
%
%      2.0000    2.5198    2.8845    3.1748    3.4200    3.6342
%
%See also
% EQ_VERTEX_DIAM, EQ_DIAM_COEFF

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

coeff = eq_vertex_diam(dim,N) .* N.^(1/dim);
%
%end function
function vertex_diam = eq_vertex_diam(dim,N)
%EQ_VERTEX_DIAM Maximum vertex diameter of EQ partition
%
%Syntax
% vertex_diam = eq_vertex_diam(dim,N);
%
%Description
% VERTEX_DIAM = EQ_VERTEX_DIAM(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) sets VERTEX_DIAM to be the maximum vertex diameter over all the regions
%    of the partition.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result VERTEX_DIAM will an array of the same size as N.
%
%Examples
%
% >> vertex_diam = eq_vertex_diam(2,10)
%
%  vertex_diam =
%
%      1.4142
%
% >> vertex_diam = eq_vertex_diam(3,1:6)
%
%  vertex_diam =
%
%       2     2     2     2     2     2
%
%See also
% EQ_DIAM_BOUND, EQ_VERTEX_DIAM_COEFF, EQ_DIAM_COEFF, EQ_REGIONS_PROPERTY

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

vertex_diam = eq_regions_property(@max_vertex_diam_of_regions,dim,N);
%
% end function
