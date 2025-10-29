function coeff = calc_dist_coeff(dim,N,min_euclidean_dist)
%CALC_DIST_COEFF Coefficient of minimum distance
%
%Syntax
% coeff = calc_dist_coeff(dim,N,min_euclidean_dist);
%
%Description
% COEFF = CALC_DIST_COEFF(dim,N,MIN_EUCLIDEAN_DIST) sets COEFF to be the
% coefficient in the expression for the lower bound on the minimum distance of
% a minimum energy point set:
%
%    MIN_EUCLIDEAN_DIST >= COEFF N^(-1/dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The argument MIN_EUCLIDEAN_DIST must be an array of real nubers of the same array size as N.
% The result COEFF will be an array of the same size as N.
%
%Notes
% The expression for the lower bound on minimum distance of a minimum r^(-s)
% energy point set on S^dim was given by [RakSZ95] for s == 0 and dim = 2,
% [Dahl78] for s == dim-1, [KuiSS04 Theorem 8] for dim-1 <= s < dim and
% [KuiS98 (1.12) p. 525] for s > dim.
%
%Examples
%
% >> N = 2:6
%
%  N =
%       2     3     4     5     6
%
% >> dist = eq_min_dist(2,N)
%
%  dist =
%
%      2.0000    1.4142    1.4142    1.4142    1.4142
%
% >> calc_dist_coeff(2,N,dist)
%
%  ans =
%
%      2.8284    2.4495    2.8284    3.1623    3.4641
%
%See also
% EQ_MIN_DIST, EQ_DIST_COEFF

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
narginchk(3,3);
%
% dim is the number of dimensions
% N is the number of regions
%
coeff =  min_euclidean_dist .* N.^(1/dim);
%
% end function
function coeff = calc_energy_coeff(dim,N,s,energy)
%CALC_ENERGY_COEFF Coefficient of second term in expansion of energy
%
%Syntax
% coeff = calc_energy_coeff(d,N,s,energy);
%
%Description
% COEFF = CALC_ENERGY_COEFF(dim,N,s,ENERGY) sets COEFF to be the coefficient of
% the second term of an expansion of ENERGY with the same form as the expansion
% of E(dim,N,s), the minimum r^(-s) energy of a set of N points on S^dim.
%
% Specifically, for s not equal to 0, COEFF is the solution to
%
% ENERGY == (SPHERE_INT_ENERGY(dim,s)/2) N^2 + COEFF N^(1+s/dim),
%
% and for s == 0 (the logarithmic potential), COEFF is the solution to
%
% ENERGY == (SPHERE_INT_ENERGY(dim,0)/2) N^2 + COEFF N LOG(N).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The argument ENERGY must an array of real numbers of the same array size as N.
% The result COEFF will be an array of the same size as N.
%
%Notes
% 1) The energy expansion is not valid for N == 1, and in particular,
%
% EQ_ENERGY_COEFF(dim,N,0,energy) := 0.
%
% 2) For s > 0, [KuiS98 (1.6) p524] has
%
% E(dim,N,s) == (SPHERE_INT_ENERGY(dim,s)/2) N^2 + COEFF N^(1+s/dim) + ...
%
% where SPHERE_INT_ENERGY(dim,s) is the energy integral of the r^(-s) potential
% on S^dim.
%
% The case s == 0 (logarithmic potential) can be split into subcases.
% For s == 0 and dim == 1, E(1,N,0) is obtained by equally spaced points on S^1,
% and the formula for the log potential for N equally spaced points on a circle
% gives
%
% E(1,N,0) == (-1/2) N LOG(N) exactly.
%
% For s == 0 and dim == 2, [SafK97 (4) p7] has
%
% E(2,N,0) == (SPHERE_INT_ENERGY(2,0)/2) N^2 + COEFF N LOG(N) + o(N LOG(N)).
%
% In general, for s == 0,
%
% E(dim,N,0) == (SPHERE_INT_ENERGY(dim,0)/2) N^2 + COEFF N LOG(N) + ...
%
% with sphere_int_energy(1,0) == 0.
%
% CALC_ENERGY_COEFF just uses this general formula for s == 0, so for s == 0 and
% dim == 1, the coefficient returned is actually the coefficient of the first
% non-zero term.
%
%Examples
%
% >> N = 2:6
%
%  N =
%
%       2     3     4     5     6
%
% >> energy = eq_energy_dist(2,N,0)
%
%  energy =
%
%     -0.6931   -1.3863   -2.7726   -4.4205   -6.2383
%
% >> coeff = calc_energy_coeff(2,N,0,energy)
%
%  coeff =
%
%     -0.2213   -0.1569   -0.2213   -0.2493   -0.2569
%
%See also
% EQ_ENERGY_DIST, EQ_ENERGY_COEFF

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
% Compute the energy coefficient: subtract the first term in the expansion of
% the minimum energy and divide by the power of N in the second term.
%
if s > 0
    %
    % The first term in the expansion of the minimum energy.
    % Ref: [KuiS98 (1.6) p524]
    %
    first_term = (sphere_int_energy(dim,s)/2) * N.^2;
    coeff = (energy-first_term) ./ (N.^(1+s/dim));
else
    %
    % Flatten N into a row vector.
    %
    shape = size(N);
    n_partitions = prod(shape);
    N = reshape(N,1,n_partitions);
    %
    % Refs for s==0, dim == 2:
    % [SafK97 (4) p. 7] [Zho95 (5.6) p. 68, 3.11 - corrected) p. 42]
    %
    first_term = (sphere_int_energy(dim,s)/2) * N.^2;
    %
    % Avoid division by zero.
    %
    coeff = zeros(size(N));
    neq1 = (N ~= 1);
    coeff(neq1) = (energy(neq1)-first_term(neq1)) ./ (N(neq1) .* log(N(neq1)));
    %
    % Reshape output to same array size as original N.
    %
    coeff = reshape(coeff,shape);
    %
end
%
% end function

function energy = sphere_int_energy(dim,s)
%SPHERE_INT_ENERGY Energy integral of r^(-s) potential
%
%Syntax
% energy = sphere_int_energy(d,s);
%
%Description
% ENERGY = SPHERE_INT_ENERGY(dim,s) sets ENERGY to be the energy integral
% on S^dim of the r^(-s) potential, defined using normalized Lebesgue measure.
%
% Ref for s > 0: [KuiS98 (1.6) p524]
% Ref for s == 0 and dim == 2: SafK97 (4) p. 7]
% For s == 0 and dim >= 2, integral was obtained using Maple:
%     energy = (1/2)*(omega(dim)/omega(dim+1)* ...
%          int(-log(2*sin(theta/2)*(sin(theta))^(dim-1),theta=0..Pi),
%     where omega(dim+1) == area_of_sphere(dim).

if s ~= 0
    energy = (gamma((dim+1)/2)*gamma(dim-s)/(gamma((dim-s+1)/2)*gamma(dim-s/2)));
elseif dim ~= 1
    energy = (psi(dim)-psi(dim/2)-log(4))/2;
else
    energy = 0;
end
%
% end function
function density = calc_packing_density(dim,N,min_euclidean_dist)
%CALC_PACKING_DENSITY Density of packing given by minimum distance
%
%Syntax
% density = calc_packing_density(dim,N,min_euclidean_dist);
%
%Description
% DENSITY = CALC_PACKING_DENSITY(dim,N,MIN_EUCLIDEAN_DIST) sets DENSITY to
% be the density of a packing of S^dim by N equal spherical caps with centers
% having minimum Euclidean distance MIN_EUCLIDEAN_DIST.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The argument MIN_EUCLIDEAN_DIST must be an array of real nubers of the same array size as N.
% The result DENSITY will be an array of the same size as N.
%
%Notes
% The packing density is defined to be the sum of the areas of the spherical
% caps of the packing, divided by the area of the unit sphere S^dim.
%
% The spherical radius of the caps in the packing is half the minimum spherical
% distance between center points. The spherical radius for N == 1 is a special
% case. It is defined to be pi.
%
%Examples
%
% >> N = 2:6
%
%  N =
%       2     3     4     5     6
%
% >> dist = eq_min_dist(2,N)
%
%  dist =
%
%      2.0000    1.4142    1.4142    1.4142    1.4142
%
% >> density = calc_packing_density(2,N,dist)
%
%  density =
%      1.0000    0.4393    0.5858    0.7322    0.8787
%
%See also
% EQ_MIN_DIST, AREA_OF_CAP, AREA_OF_SPHERE, EQ_PACKING_DENSITY

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.10 $ $Date 2005-05-27 $
% Function e2s changed name to euc2sph_dist
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
narginchk(3,3);
%
s_cap = euc2sph_dist(min_euclidean_dist)/2;
%
% The spherical radius for N == 1 is a special case. It is pi.
%
s_cap(N == 1) = pi;
density = N.*area_of_cap(dim, s_cap)./area_of_sphere(dim);
%
% end function
function coeff = eq_dist_coeff(dim,N,varargin)
%EQ_DIST_COEFF Coefficient of minimum distance of an EQ point set
%
%Syntax
% coeff = eq_dist_coeff(dim,N,options);
%
%Description
% COEFF = EQ_DIST_COEFF(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region,
% 3) finds the minimum Euclidean distance between points of the EQ point set,
% 4) sets COEFF to be the coefficient in the expression for the lower bound on
%    the minimum distance of a minimum energy point set:
%
%    DIST >= COEFF N^(-1/dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result COEFF will be an array of the same size as N.
%
% COEFF = EQ_DIST_COEFF(dim,N,'offset','extra'), for dim == 2 or dim == 3, uses
% experimental extra rotation offsets to try to maximize the minimum distance.
% For dim > 3, extra offsets are not used.
%
%Notes
% The expression for the lower bound on minimum distance of a minimum r^(-s)
% energy point set on S^dim was given by [RakSZ95] for s == 0 and dim = 2,
% [Dahl78] for s == dim-1, [KuiSS04 Theorem 8] for dim-1 <= s < dim and
% [KuiS98 (1.12) p. 525] for s > dim.
%
% Ideally eq_dist_coeff(dim,N) should tend to area_of_sphere(dim)^(1/dim) as
% N goes to infinity.
%
%Examples
%
% >> coeff = eq_dist_coeff(2,10)
%
%  coeff =
%
%      3.3250
%
% >> coeff = eq_dist_coeff(3,1:6)
%
%  coeff =
%
%      2.0000    2.5198    2.0396    2.2449    2.4183    2.5698
%
%See also
% PARTITION_OPTIONS, EQ_MIN_DIST

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
narginchk(2,4);
%
% dim is the number of dimensions
% N is the number of regions
%
dist = eq_min_dist(dim,N,varargin{:});
coeff =  dist .* N.^(1/dim);
%
% end function
function coeff = eq_energy_coeff(dim,N,s,varargin)
%EQ_ENERGY_COEFF Coefficient in expansion of energy of an EQ point set
%
%Syntax
% coeff = eq_energy_coeff(dim,N,s,options);
%
%Description
% COEFF = EQ_ENERGY_COEFF(dim,N,s)  does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region,
% 3) finds the r^(-s) energy of the EQ point set, and
% 4) sets COEFF to be the coefficient of the second term of the expansion of
%    the energy, having the same form as the expansion of E(dim,N,s),
%    the minimum r^(-s) energy of a set of N points on S^dim.
%
% Specifically, for s not equal to 0, COEFF is the solution to
%
% energy == (SPHERE_INT_ENERGY(dim,s)/2) N^2 + COEFF N^(1+s/dim),
%
% and for s == 0 (the logarithmic potential), COEFF is the solution to
%
% energy == (SPHERE_INT_ENERGY(dim,0)/2) N^2 + COEFF N LOG(N).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result COEFF will be an array of the same size as N.
%
% COEFF = EQ_ENERGY_COEFF(dim,N) uses the default value dim-1 for s.
%
% COEFF = EQ_ENERGY_COEFF(dim,N,s,'offset','extra') uses experimental extra offsets
% for S^2 and S^3 to try to minimize energy. For dim > 3, extra offsets are
% not used.
%
%Notes
% 1) The energy expansion is not valid for N == 1, and in particular,
%
% EQ_ENERGY_COEFF(dim,N,0) := 0.
%
% 2) For details of calculation of the energy coefficient,
% see HELP CALC_ENERGY_COEFF
%
%Examples
%
% >> coeff = eq_energy_coeff(2,10)
%
%  coeff =
%
%     -0.5461
%
% >> coeff = eq_energy_coeff(3,1:6)
%
%  coeff =
%
%     -0.5000   -0.5512   -0.5208   -0.5457   -0.5472   -0.5679
%
% >> coeff = eq_energy_coeff(2,1:6,0)
%
%  coeff =
%
%           0   -0.2213   -0.1569   -0.2213   -0.2493   -0.2569
%
%See also
% PARTITION_OPTIONS, EQ_ENERGY_DIST, CALC_ENERGY_COEFF

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
narginchk(2,5);
%
% dim is the number of dimensions
% N is the number of regions
%
% The default value of s is dim-1.
%
if nargin < 3
    s = dim-1;
end
energy = eq_energy_dist(dim,N,s,varargin{:});
coeff = calc_energy_coeff(dim,N,s,energy);
%
% end function

function [energy,dist] = eq_energy_dist(dim,N,s,varargin)
%EQ_ENERGY_DIST Energy and minimum distance of an EQ point set
%
% Syntax
% [energy,dist] = eq_energy_dist(dim,N,s,options);
%
%Description
% [ENERGY,DIST] = EQ_ENERGY_DIST(dim,N,s) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region,
% 3) sets ENERGY to be the r^(-s) energy of the EQ point set, and
% 4) optionally, sets DIST to be the minimum Euclidean distance between
%    points of the EQ point set.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The results ENERGY and DIST will be arrays of the same size as N.
% The result DIST is optional.
%
% [ENERGY,DIST] = EQ_ENERGY_DIST(dim,N) uses the default value dim-1 for s.
%
% [ENERGY,DIST] = EQ_ENERGY_DIST(dim,N,s,'offset','extra') uses experimental
% extra offsets for S^2 and S^3 to try to minimize energy.
% For dim > 3, extra offsets are not used.
%
%Examples
%
% >> energy = eq_energy_dist(2,10)
%
%  energy =
%     32.7312
%
% >> [energy,dist] = eq_energy_dist(3,1:6,0)
%
%  energy =
%           0   -0.6931   -1.3863   -2.7726   -4.1589   -6.2383
%
%  dist =
%      2.0000    2.0000    1.4142    1.4142    1.4142    1.4142
%
% >> [energy,dist] = eq_energy_dist(3,100,1,'offset','extra')
%
%  energy =
%     4.0042e+03
%
%  dist =
%      0.6545
%
%See also
% EQ_POINT_SET, PARTITION_OPTIONS, POINT_SET_ENERGY_DIST, EUCLIDEAN_DIST,
% MIN_DIST

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
% Check number of arguments.
%
narginchk(2,5);
%
if nargin < 3
    %
    % The default value of s is dim-1.
    %
    s = dim-1;
elseif ischar(s)
    %
    % Ensure that the user has not omitted argument s.
    %
    error('Argument s must be numeric.');
end
%
% If the option 'offset' is 'extra', then use experimental extra offsets
% for S^2 and S^3 regions to try to minimize energy.
% The default is not to use extra offsets.
%
pdefault.extra_offset =  false;
if nargin < 4
    extra_offset = pdefault.extra_offset;
else
    popt = partition_options(pdefault, varargin{:});
    extra_offset = popt.extra_offset;
end
%
% Flatten N into a row vector.
%
shape = size(N);
n_partitions = prod(shape);
N = reshape(N,1,n_partitions);

energy = zeros(size(N));
if nargout > 1
    dist = zeros(size(N));
end
for partition_n = 1:n_partitions
    points = eq_point_set(dim,N(partition_n),extra_offset);
    %
    if nargout > 1
        [energy(partition_n),dist(partition_n)] = point_set_energy_dist(points,s);
    else
        energy(partition_n) = point_set_energy_dist(points,s);
    end
end
%
% Reshape output to same array size as original N.
%
energy = reshape(energy,shape);
if nargout > 1
    dist = reshape(dist,shape);
end
%
% end function
function dist = eq_min_dist(dim,N,varargin)
%EQ_MIN_DIST Minimum distance between center points of an EQ partition
%
%Syntax
% dist = eq_min_dist(dim,N,options);
%
%Description
% DIST = EQ_MIN_DIST(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region, and
% 3) sets DIST to be the minimum Euclidean distance between points of
%    the EQ point set.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result DIST will be an array of the same size as N.
%
% DIST = EQ_MIN_DIST(dim,N,'offset','extra'), for dim == 2 or dim == 3,
% uses exerimental extra rotation offsets to try to maximize the minimum
% distance. For dim > 3, extra offsets are not used.
%
%Examples
%
% >> dist = eq_min_dist(2,10)
%
%  dist =
%      1.0515
%
% >> dist = eq_min_dist(3,1:6)
%
%  dist =
%      2.0000    2.0000    1.4142    1.4142    1.4142    1.4142
%
%See also
% PARTITION_OPTIONS, EUCLIDEAN_DIST, EQ_ENERGY_DIST

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
narginchk(2,4);
dist = eq_point_set_property(@point_set_min_dist,dim,N,varargin{:});
%
% end function
function density = eq_packing_density(dim,N,varargin)
%EQ_PACKING_DENSITY  Density of packing given by minimum distance of EQ point set
%
%Syntax
% density = eq_packing_density(dim,N,options);
%
%Description
% DENSITY = EQ_PACKING_DENSITY(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region,
% 3) finds the minimum Euclidean distance between points of the EQ point set,
% 4) sets DENSITY to be the maximum density of a packing of S^dim by equal
%    spherical caps with centers at the EQ point set.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result DENSITY will be an array of the same size as N.
%
% DENSITY = EQ_PACKING_DENSITY(dim,N,'offset','extra'), for dim == 2 or dim == 3,
% uses experimental extra rotation offsets to try to maximize the minimum
% distance. For dim > 3, extra offsets are not used.
%
%Notes
% The packing density is defined to be the sum of the areas of the spherical
% caps of the packing, divided by the area of the unit sphere S^dim.
%
% The spherical radius of the caps in the packing is half the minimum spherical
% distance between center points. The spherical radius for N == 1 is a special
% case. It is defined to be pi.
%
%Examples
%
% >> density = eq_packing_density(2,10)
%
%  density =
%
%      0.7467
%
% >> density = eq_packing_density(3,1:6)
%
%  density =
%
%      1.0000    1.0000    0.2725    0.3634    0.4542    0.5451
%
%See also
% EQ_MIN_DIST, AREA_OF_CAP, AREA_OF_SPHERE, PARTITION_OPTIONS

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
narginchk(2,4);
%
min_euclidean_dist = eq_min_dist(dim,N,varargin{:});
density = calc_packing_density(dim,N,min_euclidean_dist);
%
% end function
function property = eq_point_set_property(fhandle,dim,N,varargin)
%EQ_POINT_SET_PROPERTY Property of an EQ point set
%
%Syntax
% property = eq_point_set_property(fhandle,dim,N,options);
%
%Description
% PROPERTY = EQ_POINT_SET_PROPERTY(FHANDLE,dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
%    partition the unit sphere S^dim into N regions,
% 2) finds the EQ point set, the set of center points of each region,
% 3) calls the function defined by FHANDLE which is expected to use the
%    EQ point set to calculate the value of the result, PROPERTY.
%
% The argument FHANDLE must be a function handle. The function specified by
% FHANDLE must take as its argument a single (dim+1 by N) array, representing
% N points of S^dim, in Cartesian coordinates, and must return a single value
% based on these points.
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result PROPERTY will be an array of the same size as N.
%
% PROPERTY = EQ_POINT_SET_PROPERTY(FHANDLE,dim,N,'offset','extra'),
% for dim == 2 or dim == 3, uses exerimental extra rotation offsets to try to
% maximize the minimum distance between points of the EQ point set.
% For dim > 3, extra offsets are not used.
%
%Examples
%
% >> dist = eq_point_set_property(@point_set_min_dist,2,10)
%
%  dist =
%
%      1.0515
%
%See also
% EQ_POINT_SET, FEVAL, PARTITION_OPTIONS

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
narginchk(3,5);
%
% dim is the number of dimensions
% N is the number of regions
%
% If the option 'offset' is 'extra', then use experimental extra offsets
% for S^2 and S^3 regions to try to minimize energy.
% The default is not to use extra offsets.
%
pdefault.extra_offset =  false;
popt = partition_options(pdefault, varargin{:});
extra_offset = popt.extra_offset;
%
% Flatten N into a row vector.
%
shape = size(N);
n_partitions = prod(shape);
N = reshape(N,1,n_partitions);

property = zeros(size(N));
for partition_n = 1:n_partitions
    points = eq_point_set(dim,N(partition_n),extra_offset);
    property(partition_n) = feval(fhandle,points);
end
%
% Reshape output to same array size as original N.
%
property = reshape(property,shape);
%
% end function
function coeff = point_set_dist_coeff(points)
%POINT_SET_DIST_COEFF Coefficient of minimum distance of a point set
%
%Syntax
% coeff = point_set_dist_coeff(points);
%
%Description
% COEFF = POINT_SET_DIST_COEFF(POINTS) does the following:
% 1) finds the minimum Euclidean distance between points of the point set
%    POINTS, which should be a subset of the unit sphere S^dim, and
% 2) sets COEFF to be the coefficient in the expression for the lower bound on
%    the minimum distance of a minimum energy point set:
%
%    MIN_EUCLIDEAN_DIST >= COEFF N^(-1/dim),
%
%    where N is the number of points in POINTS.
%
% The argument POINTS must be an array of real numbers of size (dim+1 by N),
% where dim and N are positive integers.
% Each column of POINTS represents a point in R^(dim+1).
% It is assumed that point set POINTS is a subset of the unit sphere S^dim,
% but this is not checked.
%
%Notes
% Fore more details on the calculation of the coefficient of the minimum
% distance, see HELP CALC_DIST-COEFF.
%
%Examples
%
% >> x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
%
% x =
%
%      0     0     0     0
%      0     1    -1     0
%      1     0     0    -1
%
% >> coeff = point_set_dist_coeff(x)
%
% coeff =
%
%     2.8284
%
%See also
% POINT_SET_MIN_DIST, CALC_DIST_COEFF, EQ_DIST_COEFF, EQ_MIN_DIST

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-12 $
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
narginchk(1,1);
%
% dim is the dimension of S^dim as a manifold.
%
dim = size(points,1)-1;
%
% N is the number of points in the point set.
%
N = size(points,2);
%
min_euclidean_dist = point_set_min_dist(points);
coeff = calc_dist_coeff(dim,N,min_euclidean_dist);
%
% end function
function coeff = point_set_energy_coeff(points,s)
%POINT_SET_ENERGY_COEFF Coefficient in expansion of energy of a point set
%
%Syntax
% coeff = point_set_energy_coeff(points,s);
%
%Description
% COEFF = POINT_SET_ENERGY_COEFF(POINTS,s)  does the following:
% 1) finds the r^(-s) energy of the point set POINTS, and
% 2) sets COEFF to be the coefficient of the second term of the expansion of
%    the energy, having the same form as the expansion of E(dim,N,s),
%    the minimum r^(-s) energy of a set of N points on S^dim.
%
% The argument POINTS must be an array of real numbers of size (dim+1 by N),
% where dim and N are positive integers.
% Each column of POINTS represents a point in R^(dim+1).
%
% Specifically, for s not equal to 0, COEFF is the solution to
%
% energy == (SPHERE_INT_ENERGY(dim,s)/2) N^2 + COEFF N^(1+s/dim),
%
% and for s == 0 (the logarithmic potential), COEFF is the solution to
%
% energy == (SPHERE_INT_ENERGY(dim,0)/2) N^2 + COEFF N LOG(N).
%
% COEFF = POINT_SET_ENERGY_COEFF(POINTS) uses the default value dim-1 for s.
%
%Notes
% The value dim is the dimension of S^dim as a manifold. The point set POINTS
% is assumed to be a subset of R^(dim+1) but is not assumed to be a subset of
% S^dim.
%
% For details of the calculation of the energy coefficient,
% see HELP CALC_ENERGY_COEFF.
%
%Examples
%
% >> x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
%
% x =
%
%      0     0     0     0
%      0     1    -1     0
%      1     0     0    -1
%
% >> coeff = point_set_energy_coeff(x)
%
%  coeff =
%     -0.5214
%
%See also
% POINT_SET_ENERGY_DIST, CALC_ENERGY_COEFF, EQ_ENERGY_COEFF, EQ_ENERGY_DIST

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-12 $
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
narginchk(1,2);
%
% dim is the dimension of S^dim as a manifold.
%
dim = size(points,1)-1;
%
% N is the number of points in the point set.
%
N = size(points,2);
%
% The default value of s is dim-1.
%
if nargin < 2
    s = dim-1;
end

energy = point_set_energy_dist(points,s);
coeff = calc_energy_coeff(dim,N,s,energy);
%
% end function

function [energy,min_dist] = point_set_energy_dist(points,s)
%POINT_SET_ENERGY_DIST Energy and minimum distance of a point set
%
%Syntax
% [energy,min_dist] = point_set_energy_dist(points,s);
%
%Description
% [ENERGY,MIN_DIST] = POINT_SET_ENERGY_DIST(POINTS,s) sets ENERGY to be the
% energy of the r^(-s) potential on the point set POINTS, and sets MIN_DIST
% to be the minimum Euclidean distance between points of POINTS.
%
% POINTS must be an array of real numbers of size (M by N), where M and N
% are positive integers, with each of the N columns representing a point of
% R^M in Cartesian coordinates.
% The result MIN_DIST is optional.
%
% [ENERGY,MIN_DIST] = POINT_SET_ENERGY_DIST(POINTS) uses the default value of
% dim-1 for s.
%
%Notes
% The value of ENERGY for a single point is 0.
% Since this function is usually meant to be used for points on a unit sphere,
% the value of MIN_DIST for a single point is defined to be 2.
%
%Examples
%
% >> x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
%
% x =
%
%      0     0     0     0
%      0     1    -1     0
%      1     0     0    -1
%
% >> [energy,min_dist] = point_set_energy_dist(x)
%
% energy =
%
%     2.5000
%
% min_dist =
%
%     1.4142
%
%See also
% EUCLIDEAN_DIST, EQ_ENERGY_DIST, POINT_SET_MIN_DIST

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-12 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-22 $
% Documentation files renamed
% Fix nasty but obvious bug by using separate variables dist and min_dist
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(points,1)-1;
%
% The default value of s is dim-1.
%
if nargin < 2
    s = dim-1;
end
energy = 0;
if nargout > 1
    min_dist = 2;
end
n_points = size(points,2);
for i = 1:(n_points-1)
    j = (i+1):n_points;
    dist = euclidean_dist(points(:,i)*ones(1,n_points-i),points(:,j));
    energy = energy + sum(potential(s,dist));
    if nargout > 1
        min_dist = min(min_dist,min(dist));
    end
end
%
% end function

function pot = potential(s,r)
%POTENTIAL r^(-s) potential at Euclidean radius r.
%
%Syntax
% pot = potential(s,r);

switch s
    case 0
        pot = -log(r);
    otherwise
        pot = r.^(-s);
end
%
% end function
function min_dist = point_set_min_dist(points)
%POINT_SET_MIN_DIST Minimum distance between points of a point set
%
%Syntax
% min_dist = point_set_min_dist(points);
%
%Description
% MIN_DIST = POINT_SET_MIN_DIST(POINTS) sets MIN_DIST to be the minimum
% Euclidean distance between points of the point set POINTS.
%
% POINTS must be an array of real numbers of size (M by N), where M and N
% are positive integers, with each of the N columns representing a point of
% R^M in Cartesian coordinates.
%
%Notes
% Since this function is usually meant to be used for points on a unit sphere,
% the value of POINT_SET_MIN_DIST for a single point is defined to be 2.
%
%Examples
%
% >> x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
%
% x =
%
%      0     0     0     0
%      0     1    -1     0
%      1     0     0    -1
%
% >> min_dist = point_set_min_dist(x)
%
% min_dist =
%
%     1.4142
%
%See also
% EUCLIDEAN_DIST, MIN_EUCLIDEAN_DIST, POINT_SET_ENERGY_DIST

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-12 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

min_dist = 2;
n_points = size(points,2);
for i = 1:(n_points-1)
    j = (i+1):n_points;
    dist = euclidean_dist(points(:,i)*ones(1,n_points-i),points(:,j));
    min_dist = min(min_dist,min(dist));
end
%
% end function
function density = point_set_packing_density(points)
%POINT_SET_PACKING_DENSITY  Density of packing given by minimum distance of a point set
%
%Syntax
% density = point_set_packing_density(points);
%
%Description
% DENSITY = POINT_SET_PACKING_DENSITY(POINTS) does the following:
% 1) finds the minimum Euclidean distance between points of the point set
%    POINTS, and
% 2) sets DENSITY to be the density of a packing of S^dim by N equal
% spherical caps with this minimum distance.
%
% The argument POINTS must be an array of real numbers of size (dim+1 by N),
% where dim and N are positive integers.
% Each column of POINTS must represents a point of S^dim in Cartesian
% coordinates.
%
%Notes
% Because packing density is defined using spherical caps, it well defined only
% for points on S^dim. Therefore POINTS must represent a subset of S^dim.
%
% For more details on the calculation of packing density,
% see HELP CALC_PACKING_DENSITY.
%
%Examples
%
% >> x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
%
% x =
%
%      0     0     0     0
%      0     1    -1     0
%      1     0     0    -1
%
% >> density = point_set_packing_density(x)
%
%  density =
%
%      0.5858
%
%See also
% CALC_PACKING_DENSITY, EQ_MIN_DIST, AREA_OF_CAP, AREA_OF_SPHERE

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-12 $
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
narginchk(1,1);
%
% Check that points lie on the unit sphere.
%
tol = eps * 2^5;
radius = sqrt(sum(points.*points));
if (min(radius) < 1-tol) || (max(radius) > 1+tol)
    error('Point set must be a subset of the unit sphere');
end
%
% dim is the dimension of S^dim as a manifold.
%
dim = size(points,1)-1;
%
% N is the number of points in the point set.
%
N = size(points,2);
%
min_euclidean_dist = min(2,point_set_min_dist(points));
density = calc_packing_density(dim,N,min_euclidean_dist);
%
% end function
