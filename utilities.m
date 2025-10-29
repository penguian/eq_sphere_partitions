function area = area_of_cap(dim, s_cap)
%AREA_OF_CAP Area of spherical cap
%
%Syntax
% area = area_of_cap(dim, s_cap);
%
%Description
% AREA = AREA_OF_CAP(dim, S_CAP) sets AREA to be the area of an S^dim spherical
% cap of spherical radius S_CAP.
%
% The argument dim must be a positive integer.
% The argument S_CAP must be a real number or an array of real numbers.
% The result AREA will be an array of the same size as S_CAP.
%
%Notes
% S_CAP is assumed to be in the range [0, pi].
%
% The area is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% For dim <= 2, and for dim==3 (when pi/6 <= s_cap <= pi*5/6),
% AREA is calculated in closed form, using the analytic solution of
% the definite integral given in the reference.
% Otherwise, AREA is calculated using the Matlab function BETAINC,
% the incomplete Beta function ratio.
%
% Ref: [LeGS01 Lemma 4.1 p255].
%
%Examples
%
% >> a = area_of_cap(2,pi/2)
%
% a =
%
%     6.2832
%
% >> a = area_of_cap(3,0:pi/4:pi)
%
% a =
%
%          0    1.7932    9.8696   17.9460   19.7392
%
%See also
% SRADIUS_OF_CAP

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.02 $ $Date 2005-04-24 $ PL
% Use incomplete Beta function BETAINC for dim == 3,
% (when s_cap < pi/6 or s_cap > pi*5/6) and for all dim > 3.
% Use sin(s_cap).^2 in preference to (1-cos(s_cap))/2.
% $Revision 1.01 $ $Date 2005-03-16 $ PL
% Use incomplete Beta function BETAINC for dim > 8.
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

switch dim
case 1
    area = 2 * s_cap;
case 2
    area = 4*pi * sin(s_cap/2).^2;
case 3
     %
     % Flatten s_cap into a row vector.
     %
     shape = size(s_cap);
     n = prod(shape);
     s_cap = reshape(s_cap,1,n);
     area = zeros(size(s_cap));
     %
     % Near the poles, use the incomplete Beta function ratio.
     %
     pole = (s_cap < pi/6) | (s_cap > pi*5/6);
     area(pole) = area_of_sphere(dim) * betainc(sin(s_cap(pole)/2).^2,dim/2,dim/2);
     %
     % In the tropics, use closed solution to integral.
     %
     trop = s_cap(~pole);
     area(~pole) = (2*trop-sin(2*trop))*pi;

     area = reshape(area,shape);
otherwise
     area = area_of_sphere(dim) * betainc(sin(s_cap/2).^2,dim/2,dim/2);
end
%
% end function
function area = area_of_collar(dim, a_top, a_bot)
%AREA_OF_COLLAR Area of spherical collar
%
%Syntax
% area = area_of_collar(dim, a_top, a_bot);
%
%Description
% AREA = AREA_OF_COLLAR(dim, A_TOP, A_BOT) sets AREA to be the area of
% an S^dim spherical collar specified by A_TOP, A_BOT, where
% A_TOP is top (smaller) spherical radius,
% A_BOT is bottom (larger) spherical radius.
%
% The argument dim must be a positive integer.
% The arguments A_TOP and A_BOT must be real numbers or arrays of real numbers,
% with the same array size.
% The result AREA will be an array of the same size as A_TOP.
%
%Notes
% A_TOP and A_BOT are assumed to be in the range [0, pi].
%
% The area is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% Ref: [LeGS01 Lemma 4.1 p255].
%
%Examples
%
% >> area = area_of_collar(2,0:2,1:3)
%
% area =
%
%     2.8884    6.0095    3.6056
%
%See also
% AREA_OF_CAP, AREA_OF_SPHERE

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

area = area_of_cap(dim, a_bot) - area_of_cap(dim, a_top);
%
% end function
function area = area_of_ideal_region(dim,N)
%AREA_OF_IDEAL_REGION Area of one region of an EQ partition
%
%Syntax
% area = area_of_ideal_region(dim,N);
%
%Description
% AREA = AREA_OF_IDEAL_REGION(dim,N) sets AREA to be the area of one of N equal
% area regions on S^dim, that is 1/N times AREA_OF_SPHERE(dim).
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result AREA will be an array of the same size as N.
%
%Examples
%
% >> area = area_of_ideal_region(3,1:6)
%
% area =
%
%    19.7392    9.8696    6.5797    4.9348    3.9478    3.2899
%
%See also
% AREA_OF_SPHERE

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

area = area_of_sphere(dim)./N;
%
% end function
function area = area_of_sphere(dim)
%AREA_OF_SPHERE Area of sphere
%
%Syntax
% area = area_of_sphere(dim);
%
%Description
% AREA = AREA_OF_SPHERE(dim) sets AREA to be the area of the sphere S^dim,
%
% The argument dim must be a positive integer or an array of positive integers.
% The result AREA will be an array of the same size as dim.
%
%Notes
% The area of S^dim is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% Ref: [Mue98] p39.
%
%Examples
%
% >> area = area_of_sphere(1:7)
%
% area =
%
%     6.2832   12.5664   19.7392   26.3189   31.0063   33.0734   32.4697
%
%See also
% AREA_OF_CAP, VOLUME_OF_BALL

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

power = (dim+1)/2;
area = (2*pi.^power./gamma(power));
%
% end function
function s = cart2polar2(x)
%CART2POLAR2 Convert from Cartesian to spherical coordinates on sphere S^2
%
%Syntax
% s = cart2polar2(x);
%
%Description
% S = CART2POLAR2(X) sets S to be the spherical polar coordinates of the points
% represented by the Cartesian coordinates X:
% S = [phi;theta]: phi in [0, 2*pi), theta in [0, pi].
%
% The argument X must be an array of real numbers of size (3 by N), where N is
% any positive integer. The result S will be an array of size (2 by N).
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
% >> s = cart2polar2(x)
%
% s =
%          0    1.5708    4.7124         0
%          0    1.5708    1.5708    3.1416
%
%Note
% CART2POLAR2(X) projects any X in R^3 onto the sphere S^2 via a line through
% the origin. The origin [0 0 0]' is itself projected onto a point on the
% equator such that
%
%     POLAR2CART(CART2POLAR2([0 0 0]')) == [1 0 0]'.
%
%See also
% POLAR2CART

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-12 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.10 $ $Date 2005-05-28 $
% Use cart2sph
% Change name from x2s2 to cart2polar2
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

[phi, theta] = cart2sph(x(1,:),x(2,:),x(3,:));
s = [mod(phi, 2*pi); pi/2-theta];
%Recursive Zonal Equal Area Sphere Partitioning: Utilities
%
% Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox.
% Release 1.12 2024-10-14
%
%Functions
%=========
%
%  area_of_cap            Area of spherical cap
%  area_of_collar         Area of spherical collar
%  area_of_ideal_region   Area of one region of an EQ partition
%  area_of_sphere         Area of sphere
%  cart2polar2            Convert Cartesian to spherical polar coordinates on S^2
%  euc2sph_dist           Convert Euclidean to spherical distance
%  euclidean_dist         Euclidean distance between two points
%  fatcurve               Create a parameterized cylindrical surface
%  ideal_collar_angle     Ideal angle for spherical collars of an EQ partition
%  polar2cart             Convert spherical polar to Cartesian coordinates
%  sph2euc_dist           Convert spherical to Euclidean distance
%  spherical_dist         Spherical distance between two points on the sphere
%  sradius_of_cap         Spherical radius of spherical cap of given area
%  volume_of_ball         Volume of the unit ball

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-14 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from e2s to euc2sph_dist
% Function changed name from s2e to sph2euc_dist
% Function changed name from s2x to polar2cart
% Function changed name from x2s2 to cart2polar2
% Add new function fatcurve
% Add new function haslight
% Clean up descriptions
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

function s = euc2sph_dist(e)
%EUC2SPH_DIST Convert Euclidean to spherical distance
%
%Syntax
% s = e2s(e);
%
%Description
% S = EUC2SPH_DIST(E) converts the Euclidean distance E to the spherical
% distance S, using a formula which is valid for the unit sphere in all
% dimensions.
%
% The argument E must be a real number or an array of real numbers.
% The result S will be an array of the same size as E.
%
%Note
% The argument E is assumed to satsify abs(E) <= 2.
%
%Examples
%
% >> s = euc2sph_dist(2)
%
%  s =
%      3.1416
%
% >> s = euc2sph_dist(0:0.5:2)
%
%  s =
%           0    0.5054    1.0472    1.6961    3.1416
%
% >> s = euc2sph_dist(-2)
%
%  s =
%     -3.1416
%
%See also
% SPH2EUC_DIST, EUCLIDEAN_DIST, SPHERICAL_DIST

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Change name from e2s to euc2sph_dist
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

s = 2*asin(e/2);
%
% end function
function distance = euclidean_dist(x,y)
%EUCLIDEAN_DIST Euclidean distance between two points in Cartesian coordinates
%
%Syntax
% distance = euclidean_dist(x,y);
%
%Description
% DISTANCE = EUCLIDEAN_DIST(X,Y) sets DISTANCE to be the Euclidean distance
% between the two points X and Y.
%
% The arguments X and Y must be arrays of the same size, M by N, where M and N
% are positive integers. Each of X and Y is assumed to represent N points in
% R^M, in Cartesian coordinates.
% The result DISTANCE will be a 1 by N array.
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
% >> y = [[0 -0.5 0.866]' [0 0.866 0.5]' [0 -0.866 -0.5]' [0 0.5 -0.866]']
%
% y =
%
%          0         0         0         0
%    -0.5000    0.8660   -0.8660    0.5000
%     0.8660    0.5000   -0.5000   -0.8660
%
% >> distance = euclidean_dist(x,y)
%
% distance =
%
%     0.5176    0.5176    0.5176    0.5176
%
%See also
% SPHERICAL_DIST, E2S, S2E

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

distance = sqrt(sum((x-y).^2));
%
% end function
function [X,Y,Z] = fatcurve(c,r)
%FATCURVE Create a parameterized cylindrical surface at radius r from curve c
%
%Syntax
% [X,Y,Z] = fatcurve(c,r);
%
%Description
% [X,Y,Z] = FATCURVE(C, R) sets X, Y and Z to be the coordinates of a
% cylindrical surface at radius R from curve C. This function is intended
% for use with the Matlab function SURF to illustrate curves in R^3.
%
%Examples
%
% >> N=5;
% >> phi = linspace(0, pi/5, N);
% >> theta = zeros(1, N);
% >> s = [theta; phi];
% >> c = polar2cart(s)
% >> r = 0.1
% >> [X,Y,Z] = fatcurve(c,r)
%
% c =
%
%          0    0.1564    0.3090    0.4540    0.5878
%          0         0         0         0         0
%     1.0000    0.9877    0.9511    0.8910    0.8090
%
% r =
%
%     0.1000
%
% X =
%
%          0    0.0055    0.0078    0.0055    0.0000   -0.0055   -0.0078   -0.0055   -0.0000
%     0.1564    0.1729    0.1798    0.1729    0.1564    0.1399    0.1331    0.1399    0.1564
%     0.3090    0.3361    0.3473    0.3361    0.3090    0.2820    0.2707    0.2820    0.3090
%     0.4540    0.4909    0.5062    0.4909    0.4540    0.4170    0.4017    0.4170    0.4540
%     0.5878    0.6247    0.6400    0.6247    0.5878    0.5508    0.5355    0.5508    0.5878
%
% Y =
%
%     0.1000    0.0707    0.0000   -0.0707   -0.1000   -0.0707   -0.0000    0.0707    0.1000
%     0.1000    0.0707    0.0000   -0.0707   -0.1000   -0.0707   -0.0000    0.0707    0.1000
%     0.1000    0.0707    0.0000   -0.0707   -0.1000   -0.0707   -0.0000    0.0707    0.1000
%     0.1000    0.0707    0.0000   -0.0707   -0.1000   -0.0707   -0.0000    0.0707    0.1000
%     0.1000    0.0707    0.0000   -0.0707   -0.1000   -0.0707   -0.0000    0.0707    0.1000
%
% Z =
%
%     1.0000    1.0705    1.0997    1.0705    1.0000    0.9295    0.9003    0.9295    1.0000
%     0.9877    1.0564    1.0849    1.0564    0.9877    0.9189    0.8905    0.9189    0.9877
%     0.9511    1.0164    1.0434    1.0164    0.9511    0.8857    0.8587    0.8857    0.9511
%     0.8910    0.9513    0.9763    0.9513    0.8910    0.8307    0.8057    0.8307    0.8910
%     0.8090    0.8693    0.8943    0.8693    0.8090    0.7487    0.7238    0.7487    0.8090

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.1 $ $Date 2025-08-09 $
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-14 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Flesh out description and examples
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(c,1);
if dim ~= 3
    error('Function fatcurve called with dim == %d but dim must be 3.', dim)
end
n = size(c,2);
m = 8;
h = 0:1/m:1;
phi = h*2*pi;
X = zeros(n,m+1);
Y = X;
Z = Y;
for k = 1:n-1
    u = c(:,k+1)-c(:,k);
    M = null(u');
    if size(M,2) ~= 2
        fprintf('size(M,2) == %d\n',size(M,2));
        disp(M);
        disp(c);
	return;
    end
    v = M(:,1);
    w = cross(u,v);
    w = w/norm(w);
    if k > 1
        minindex = 0;
        mindist = 2;
        for j = 1:m
            dist = norm(v - circ(:,j));
            if dist < mindist
                mindist = dist;
                minindex = j;
            end
        end
        offs = phi(minindex);
        circ = v*cos(phi-offs) + w*sin(phi-offs);
    else
        circ = v*cos(phi) + w*sin(phi);
    end
    XYZ = c(:,k)*ones(size(phi)) + r*circ;
    X(k,:) = XYZ(1,:);
    Y(k,:) = XYZ(2,:);
    Z(k,:) = XYZ(3,:);
end
XYZ = c(:,n)*ones(size(phi)) + r*circ;
X(n,:) = XYZ(1,:);
Y(n,:) = XYZ(2,:);
Z(n,:) = XYZ(3,:);
function angle = ideal_collar_angle(dim,N)
%IDEAL_COLLAR_ANGLE The ideal angle for spherical collars of an EQ partition
%
%Syntax
% angle = ideal_collar_angle(dim,N);
%
%Description
% ANGLE = IDEAL_COLLAR_ANGLE(dim,N) sets ANGLE to the ideal angle for the
% spherical collars of an EQ partition of the unit sphere S^dim into N regions.
%
% The argument dim must be a positive integer.
% The argument N must be a positive integer or an array of positive integers.
% The result ANGLE will be an array of the same size as N.
%
%Notes
% The ideal collar angle is determined by the side of a dim-dimensional
% hypercube of the same volume as the area of a single region of an N region
% equal area partition of S^dim.
%
% Since the EQ partition for N < 3 has no spherical collars,
% the recursive zonal equal area sphere partitioning algorithm does not use
% ideal_collar_angle(dim,N) for N < 3.
%
%Examples
%
% >> angle = ideal_collar_angle(2,10)
%
%  angle =
%
%      1.1210
%
% >> angle = ideal_collar_angle(3,1:6)
%
%  angle =
%
%      2.7026    2.1450    1.8739    1.7025    1.5805    1.4873
%
%See also
% AREA_OF_IDEAL_REGION

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

angle = area_of_ideal_region(dim,N).^(1/dim);
%
% end function
function x = polar2cart(s)
%POLAR2CART Convert spherical polar to Cartesian coordinates
%
%Syntax
% x = polar2cart(s);
%
%Description
% X = POLAR2CART(S) sets X to be the Cartesian coordinates of the points
% represented by the spherical polar coordinates S.
%
% S is assumed to be an array of size (dim by N) representing N points of
% S^dim in spherical polar coordinates, where dim and N are positive integers.
% N will be an array of size (dim+1 by N).
%
%Examples
%
% >> s = [[0 0]' [pi/2 pi/2]' [3*pi/2 pi/2]' [0 pi]']
%
% s =
%          0    1.5708    4.7124         0
%          0    1.5708    1.5708    3.1416
%
% >> x = polar2cart(s)
%
% x =
%          0    0.0000   -0.0000    0.0000
%          0    1.0000   -1.0000         0
%     1.0000    0.0000    0.0000   -1.0000
%
%See also
% CART2POLAR2

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Change name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(s,1);
n=size(s,2);
x = zeros(dim+1,n);
sinprod  = 1;
for k = dim:-1:2
    x(k+1,:) = sinprod.*cos(s(k,:));
    sinprod  = sinprod.*sin(s(k,:));
end
x(2,:)=sinprod.*sin(s(1,:));
x(1,:)=sinprod.*cos(s(1,:));
r = sqrt(sum(x.^2));
mask = (r ~= 1);
if  size(r(mask),2) > 0
    x(:,mask) = x(:,mask)./(ones(dim+1,1)*r(mask));
end
%
% end function
function e = sph2euc_dist(s)
%SPHE2EUC_DIST Convert spherical to Euclidean distance
%
%Syntax
% e = sph2euc_dist(s);
%
%Description
% E = SPHE2EUC_DIST (S) converts the spherical distance S to Euclidean
% distance E, using a formula which is valid for the unit sphere in all
% dimensions.
%
% The argument S must be a real number or an array of real numbers.
% The result E will be an array of the same size as S.
%
%Note
% The argument S is assumed to satsify abs(S) <= pi.
%
%Examples
%
% >> e = sph2euc_dist(pi)
%
% e =
%
%      2
%
% >> e = sph2euc_dist(0:pi/4:pi)
%
% e =
%
%          0    0.7654    1.4142    1.8478    2.0000
%
%See also
% EUC2SPH_DIST, EUCLIDEAN_DIST, SPHERICAL_DIST

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Change name from s2e to sph2euc_dist
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

e = 2*(sin(s/2));
%
% end function
function a_dist = spherical_dist(x,y)
%SPHERICAL_DIST Spherical distance between two points on the sphere
%
%Syntax
% a_dist = spherical_dist(x,y);
%
%Description
% A_DIST = SPHERICAL_DIST(X,Y) sets A_DIST to be the spherical distance
% between the two points X and Y.
%
% The arguments X and Y must be arrays of the same size, M by N, where M and N
% are positive integers. Each of X and Y is assumed to represent N points in
% R^M, in Cartesian coordinates.
% The result A_DIST will be a 1 by N array.
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
% >> y = [[0 -0.5 0.866]' [0 0.866 0.5]' [0 -0.866 -0.5]' [0 0.5 -0.866]']
%
% y =
%
%          0         0         0         0
%    -0.5000    0.8660   -0.8660    0.5000
%     0.8660    0.5000   -0.5000   -0.8660
%
% >> a_dist = spherical_dist(x,y)
%
% a_dist =
%
%     0.5236    0.5236    0.5236    0.5236
%
%See also
% EUCLIDEAN_DIST, E2S, S2E

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

a_dist = acos(sum(x.*y));
%
% end function
function s_cap = sradius_of_cap(dim, area)
%SRADIUS_OF_CAP Spherical radius of spherical cap of given area
%
%Syntax
% s_cap = sradius_of_cap(dim, area);
%
%Description
% S_CAP = SRADIUS_OF_CAP(dim, AREA) sets S_CAP to be the spherical radius of
% an S^dim spherical cap of area AREA.
%
% The argument dim must be a positive integer.
% The argument AREA must be a real number or an array of real numbers.
% The result S_CAP will be an array of the same size as AREA.
%
%Notes
% S_CAP is assumed to be in the range [0, pi].
%
% The area is defined via the Lebesgue measure on S^dim inherited from
% its embedding in R^(dim+1).
%
% For dim <= 2, S_CAP is calculated in closed form.
% Otherwise, S_CAP is approximated using the Matlab function FZERO.
%
% Ref: [LeGS01 Lemma 4.1 p255].
%
%Examples
%
% >> s_cap = sradius_of_cap(2,area_of_sphere(2)/2)
%
% s_cap =
%     1.5708
%
% >> s_cap = sradius_of_cap(3,(0:4)*area_of_sphere(3)/4)
%
% s_cap =
%          0    1.1549    1.5708    1.9867    3.1416
%
%See also
% FZERO, AREA_OF_CAP

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.1 $ $Date 2025-08-09 $
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.02 $ $Date 2005-04-25 $
% Use asin rather than acos to avoid subtraction.
% Use symmetry to avoid loss of accuracy in the Southern hemisphere.
% Remove check for when area is close to area of sphere.
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

switch dim
case 1
    s_cap = area/2;
case 2
    s_cap = 2*asin(sqrt(area/pi)/2);
otherwise
    %
    % Flatten area into a row vector.
    %
    shape = size(area);
    n = prod(shape);
    area = reshape(area,1,n);
    s_cap = zeros(size(area));
    for k = 1:n
        ak = area(k);
        as = area_of_sphere(dim);
        if ak >= as
            s_cap(k) = pi;
        else
            if (2*ak > as)
                ak = as - ak;
                flipped = true;
            else
                flipped = false;
            end
            area_diff = @(s) area_of_cap(dim, s) - ak;
            sk = fzero(area_diff,[0,pi]);

            if flipped
                s_cap(k) = pi - sk;
            else
                s_cap(k) = sk;
            end
        end
    end
    %
    % Reshape output to same array size as original area.
    %
    s_cap = reshape(s_cap,shape);
end
% end function
function volume = volume_of_ball(dim)
%VOLUME_OF_BALL Volume of the unit ball
%
%Syntax
% volume = volume_of_ball(dim);
%
%Description
% VOLUME = VOLUME_OF_BALL(dim) sets VOLUME to be the volume of the unit ball
% B^dim in R^dim which is enclosed by the sphere S^(dim-1).
%
% The argument dim must be a positive integer or an array of positive integers.
% The result VOLUME will be an array of the same size as dim.
%
%Notes
% The volume of B^dim is defined via the Lebesgue measure on R^dim.
%
% Ref: [WeiMW].
%
%Examples
%
% >> volume = volume_of_ball(1:7)
%
% volume =
%
%     2.0000    3.1416    4.1888    4.9348    5.2638    5.1677    4.7248
%
%See also
% AREA_OF_SPHERE

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

volume = area_of_sphere(dim-1)./dim;
%
% end function
