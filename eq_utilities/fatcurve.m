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
    fprintf('Error: dim == %d but should be 3.\n', dim)
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
