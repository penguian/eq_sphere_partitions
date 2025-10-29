function expanded_region = expand_region_for_diam(region)
%EXPAND_REGION_FOR_DIAM The set of 2^d vertices of a region
%
% Expand a region from the 2 vertex definition to the set of 2^dim vertices
% of the pseudo-region of a region, so that the Euclidean diameter of a region
% is approximated by the diameter of this set.
%
%Syntax
% expanded_region = expand_region_for_diam(region);
%
%Examples
%
% >> regions = eq_regions(2,10);
% >> region = regions(:,:,2)
%
% region =
%
%          0    1.5708
%     0.6435    1.5708
%
% >> expanded_region = expand_region_for_diam(region)
%
% expanded_region =
%
%          0    1.5708         0    1.5708
%     0.6435    0.6435    1.5708    1.5708
%
%See also
% PSEUDO_REGION_FOR_DIAM, MAX_VERTEX_DIAM_OF_REGIONS, EQ_VERTEX_DIAM

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2025-10-29 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(region,1);
if dim > 1
    s_top = region(dim,1);
    s_bot = region(dim,2);
    region_1 = expand_region_for_diam(region(1:dim-1,:));
    expanded_region = [append(region_1, s_top), append(region_1, s_bot)];
else
    expanded_region = pseudo_region_for_diam(region);
end
%
% end function

function result = append(matrix,value)
% Append a coordinate value to each column of a matrix.
%
% result = append(matrix,value);

result = [matrix; ones(1,size(matrix,2))*value];
%
% end function
function diam_bound = max_diam_bound_of_regions(regions)
%MAX_DIAM_BOUND_OF_REGIONS The maximum diameter bound in an array of regions
%
%Syntax
% diam_bound = max_diam_bound_of_regions(regions);
%
%Examples
%
% >> regions = eq_regions(2,10);
% >> max_diam_bound_of_regions(regions)
%
% ans =
%
%     1.6733
%
%See also
% MAX_VERTEX_DIAM_OF_REGIONS, EQ_DIAM_BOUND

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2025-10-29 $
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-04-28 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function s2e changed name to sph2euc_dist
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(regions,1);
if dim == 1
    diam_bound = diam_bound_region(regions(:,:,1));
else
    colatitude = -inf*ones(dim-1,1);
    diam_bound = 0;
    N = size(regions,3);
    for region_n = 1:N
        top = regions(:,1,region_n);
        if norm(top(2:dim)-colatitude) ~= 0
            colatitude = top(2:dim);
            diam_bound = max(diam_bound,diam_bound_region(regions(:,:,region_n)));
        end
    end
end
%
% end function

function diam_bound = diam_bound_region(region)
% Calculate the per-region bound on the Euclidean diameter of a region.
%
% diam_bound = diam_bound_region(region)

tol = eps*2^5;
pseudo_region = pseudo_region_for_diam(region);
dim = size(pseudo_region,1);
top = pseudo_region(:,1);
bot = pseudo_region(:,2);
s = bot(dim)-top(dim);
e = sph2euc_dist(s);
if dim == 1
    diam_bound = e;
else
    max_sin = max(sin(top(dim)),sin(bot(dim)));
    if (top(dim) <= pi/2) && (bot(dim) >= pi/2)
        max_sin = 1;
    end
    if (abs(top(dim)) < tol) || (abs(pi-bot(dim)) < tol)
        diam_bound = 2*max_sin;
    else
        region_1 = [top(1:dim-1),bot(1:dim-1)];
        diam_bound_1 = max_sin*diam_bound_region(region_1);
        diam_bound = min(2,sqrt(e^2+diam_bound_1^2));
    end
end
%
% end function
function vertex_diam = max_vertex_diam_of_regions(regions)
%MAX_VERTEX_DIAM_OF_REGIONS The max vertex diameter in a cell array of regions
%
%Syntax
% vertex_diam = max_vertex_diam_of_regions(regions);
%
%Examples
%
% >> regions = eq_regions(2,10);
% >> max_vertex_diam_of_regions(regions)
%
% ans =
%
%     1.4142
%
%See also
% MAX_DIAM_BOUND_OF_REGIONS, EQ_VERTEX_DIAM

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2025-10-29 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

dim = size(regions,1);
if dim == 1
    vertex_diam = vertex_diam_region(regions(:,:,1));
else
    colatitude = -inf*ones(dim-1,1);
    vertex_diam = 0;
    N = size(regions,3);
    for region_n = 1:N
        top = regions(:,1,region_n);
        if norm(top(2:dim)-colatitude) ~= 0
            colatitude = top(2:dim);
            vertex_diam = max(vertex_diam,vertex_diam_region(regions(:,:,region_n)));
        end
    end
end
%
% end function

function diam = vertex_diam_region(region)
% Calculate the Euclidean diameter of the set of 2^dim vertices
% of the pseudo-region of a region.
%
% diam = vertex_diam_region(region);

expanded_region = expand_region_for_diam(region);
dim = size(expanded_region,1);
full = size(expanded_region,2);
half = floor(full/2);
top = expanded_region(:,1);
bot = expanded_region(:,full);
diam = 0;
if sin(top(dim)) > sin(bot(dim))
    for point_n_1 = 1:2:half
        for point_n_2 = point_n_1+1:2:full
            x1 = polar2cart(expanded_region(:,point_n_1));
            x2 = polar2cart(expanded_region(:,point_n_2));
            diam = max(diam,euclidean_dist(x1,x2));
        end
    end
else
    for point_n_1 = full:-2:half+1
        for point_n_2 = point_n_1-1:-2:1
            x1 = polar2cart(expanded_region(:,point_n_1));
            x2 = polar2cart(expanded_region(:,point_n_2));
            diam = max(diam,euclidean_dist(x1,x2));
        end
    end
end
%
% end function
function pseudo_region = pseudo_region_for_diam(region)
%PSEUDO_REGION_FOR_DIAM Two points which maximize the vertex diameter of a region
%
%Syntax
% pseudo_region = pseudo_region_for_diam(region)
%
%Examples
%
% >> regions = eq_regions(2,10);
% >> region = regions(:,:,2)
%
% region =
%
%          0    1.5708
%     0.6435    1.5708
%
% >> pseudo_region = pseudo_region_for_diam(region)
%
% pseudo_region =
%
%          0    1.5708
%     0.6435    1.5708
%
%See also
% EXPAND_REGION_FOR_DIAM, MAX_VERTEX_DIAM_OF_REGIONS, EQ_VERTEX_DIAM

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2025-10-29 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

tol = eps*2^5;
phi_top = region(1,1);
phi_bot = region(1,2);
if phi_bot == 0
    phi_bot = 2*pi;
end
if (mod(phi_bot - phi_top, 2*pi) < tol) || (mod(phi_bot - phi_top, 2*pi) > pi)
    phi_bot = phi_top + pi;
end
pseudo_region = region;
pseudo_region(1,2) = phi_bot;
%
% end function
