%Recursive Zonal Equal Area Sphere Partitioning: Illustrations
%
% Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox.
% Release 1.12 2024-10-16
%
%Functions by category
%=====================
%
% Illustration of EQ Partitions of S^2 or S^3
%  show_s2_partition      3D illustration of an EQ partition of S^2
%  project_s2_partition   Use projection to illustrate an EQ partition of S^2
%  project_s3_partition   Use projection to illustrate an EQ partition of S^3
%
% Illustration of point sets on S^2 or S^3
%  show_r3_point_set      3D illustration of a point set
%  project_point_set      Use projection to illustrate a point set of
%                             S^2 or S^3
%
% Illustration options
%  illustration_options   Options for illustrations of EQ partitions
%
% Illustration utilities
%  haslight               Check if axis handle has a light attached

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
function lit = haslight(axish)
%HASLIGHT Check if axis handle has a light attached
%
%Syntax
% lit = haslight(axish);
%
%Description
% LIT = HASLIGHT(AXISH) sets LIT to true if the axis specified by axis
% handle AXISH has a child of type 'light', and sets LIT to false otherwise.
%
%Examples
%
% >> lit = haslight(gca)
%
% lit =
%
%   logical
%
%    0

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-14 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

lit = false;
ch = get(axish,'Children');

for k = 1:length(ch)
    type = get(ch(k),'Type');
    if strcmp(type,'light')
        lit = true;
        return;
    end
end
%
%end function
function gopt = illustration_options(gdefault, varargin)
%ILLUSTRATION_OPTIONS Options for illustrations of EQ partitions
%
%Syntax
% gopt = illustration_options(gdefault,options);
%
%Description
% GOPT = ILLUSTRATION_OPTIONS(GDEFAULT,options) collects illustration options,
% specified as name, value pairs, and places these into the structure GOPT.
% The structure GDEFAULT is used to define default option values.
%
% The structures gdefault and gopt may contain the following fields:
% fontsize:      numeric
% long_title:    boolean
% stereo:        boolean
% show_points:   boolean
% show_sphere:   boolean
% show_surfaces: boolean
%
% The following illustration options are available.
%
% 'fontsize':    Font size used in titles.
%      number    Assigns number to field gopt.fontsize.
%
% 'title':       Length of titles.
%     'long':    Long titles.
%                Sets gopt.show_title to true.
%                Sets gopt.long_title to true.
%     'short':   Short titles.
%                Sets gopt.show_title to true.
%                Sets gopt.long_title to false.
%     'none':    No titles.
%                Sets gopt.show_title to false.
%                Sets gopt.long_title to false.
%     'show':    Show default titles.
%                Sets gopt.show_title to true.
%     'hide':    Same as 'none'.
%
% 'proj':        Projection from the sphere to the plane R^2 or the space R^3.
%     'stereo':  Stereographic projection from the sphere to the whole space.
%                Sets gopt.stereo to true.
%     'eqarea':  Equal area projection from the sphere to the disk or ball.
%                Sets gopt.stereo to false.
%
% 'points':      Show or hide center points of regions.
%     'show':    Show center points of regions.
%                Sets gopt.show_points to true.
%     'hide':    Hide center points of regions.
%                Sets gopt.show_points to false.
%
% 'sphere':      Show or hide the sphere S^2.
%     'show':    Show sphere.
%                Sets gopt.show_sphere to true.
%     'hide':    Hide sphere.
%                Sets gopt.show_sphere to false.
%
% 'surf':        Show or hide surfaces of regions of a partition of S^3.
%     'show':    Show surfaces of regions.
%                Sets gopt.show_surfaces to true.
%     'hide':    Hide surfaces of regions.
%                Sets gopt.show_surfaces to false.
%
%Examples
%
% >> gdefault.fontsize=14;
% >> gopt = illustration_options(gdefault,'proj','stereo')
%
% gopt =
%
%     fontsize: 14
%       stereo: 1
%
% >> gopt = illustration_options(gdefault,'proj','stereo','fontsize',12)
%
% gopt =
%
%    fontsize: 12
%      stereo: 1

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

gopt = gdefault;
nargs = length(varargin);
nopts = floor(nargs/2);
opt_args = varargin(1:2:2*nopts-1);
for k=1:nopts
    if ~ischar([opt_args{k}])
        fprintf('Option names must be character strings\n');
        option_error(varargin{:});
    end
end
opt_vals = varargin(2:2:2*nopts);

option_name = 'fontsize';
pos = find(strcmp(option_name,opt_args));
if ~isempty(pos)
    if (isscalar(pos))
        gopt.fontsize = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
end

option_name = 'title';
pos = find(strcmp(option_name,opt_args));
if ~isempty(pos)
    if (isscalar(pos))
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'long'
        gopt.show_title = true;
        gopt.long_title = true;
    case 'short'
        gopt.show_title = true;
        gopt.long_title = false;
    case 'none'
        gopt.show_title = false;
        gopt.long_title = false;
    case 'hide'
        gopt.show_title = false;
        gopt.long_title = false;
    case 'show'
        gopt.show_title = true;
    otherwise
        value_error(value,varargin{:});
    end
end

option_name = 'proj';
pos = find(strcmp(option_name,opt_args));
if ~isempty(pos)
    if (isscalar(pos))
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'stereo'
        gopt.stereo = true;
    case 'eqarea'
        gopt.stereo = false;
    otherwise
        value_error(value,varargin{:});
    end
end

option_name = 'points';
pos = find(strcmp(option_name,opt_args));
if ~isempty(pos)
    if (isscalar(pos))
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'show'
        gopt.show_points = true;
    case 'hide'
        gopt.show_points = false;
    otherwise
        value_error(value,varargin{:});
    end
end

option_name = 'surf';
pos = find(strcmp(option_name,opt_args));
if ~isempty(pos)
    if (isscalar(pos))
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'show'
        gopt.show_surfaces = true;
    case 'hide'
        gopt.show_surfaces = false;
    otherwise
        value_error(value,varargin{:});
    end
end

option_name = 'sphere';
pos = find(strcmp(option_name,opt_args));
if ~isempty(pos)
    if (isscalar(pos))
        value = opt_vals{pos};
    else
        duplicate_error(option_name,varargin{:});
    end
    switch value
    case 'show'
        gopt.show_sphere = true;
    case 'hide'
        gopt.show_sphere = false;
    otherwise
        value_error(value,varargin{:});
    end
end
%
% end function

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
error('Please check "help illustration_options" for options');
%
% end function
function project_point_set(points,varargin)
%PROJECT_POINT_SET Use projection to illustrate a point set of S^2 or S^3
%
%Syntax
% project_point_set(points,options);
%
%Description
% PROJECT_POINT_SET(POINTS,OPTIONS) uses projection to illustrate a point set
% of S^2 or S^3, represented by the Cartesian coordinates POINTS.
%
% The argument POINTS must be an array of real numbers of size (3 by N) or
% (4 by N), where N is a positive integer.
%
% PROJECT_POINT_SET(N,options) also recognizes a number of illustration
% options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% PROJECT_POINT_SET(N,'fontsize',size)
% Font size used in titles (numeric, default 16).
%
% PROJECT_POINT_SET(N,'title','show')
% PROJECT_POINT_SET(N,'title','hide')
% Show or hide titles (default 'hide').
%
% PROJECT_POINT_SET(N,'proj','stereo')
% PROJECT_POINT_SET(N,'proj','eqarea')
% Use stereographic or equal area projection (default 'stereo').
%
%Notes
% The points are assumed to all lie on the unit sphere S^dim, where dim == 2 or
% dim == 3. The first point POINTS(:,1) should be the North pole [1,0,0]' or
% [1,0,0,0]'.
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
% >> project_point_set(x)
%
%See also
% ILLUSTRATION_OPTIONS

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

gdefault.fontsize = 16;
gdefault.show_title  = false;
gdefault.stereo = true;

gopt = illustration_options(gdefault, varargin{:});

dim = size(points,1)-1;
if dim ~= 2 && dim ~= 3
    error('project_point_set(points_x): points_x must be a point set in R^3 or or R^4');
end
N = size(points,2);

if gopt.stereo
    projection = @x2stereo;
else
    projection = @x2eqarea;
end

switch dim
case 2
    t = projection(points);
    r = real(pi-acos(points(end,N-size(t,2)+1:end)));
    limit = pi;
    if N <= 4
        c = 'k';
    else
        c = r;
    end
    s = ceil(40*r/limit)+8;
    h = scatter(t(1,:),t(2,:),s,c,'filled');
    if N < 256
        set(h,'MarkerEdgeColor',[0.5,0.5,0.5]);
    else
        set(h,'MarkerEdgeColor','none');
    end
    axis equal
    grid off
    axis off
case 3
    t = projection(points);
    r = real(pi-acos(points(end,N-size(t,2)+1:end)));
    if gopt.stereo
        limit = pi;
    else
        limit = 2*pi;
    end
    s = (r+1)/(limit*10);
    [X,Y,Z] = sphere;

    for k = 1:size(t,2)
        surf(t(1,k)+s(k)*X,t(2,k)+s(k)*Y,t(3,k)+s(k)*Z,r(k)*ones(size(X)),...
             'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
        axis equal
        hold on
    end
    grid off
    axis off
    camlight right
end

if gopt.show_title
    titlestr = sprintf(...
        'Point set containing %d points.',N);
    title(titlestr,'FontWeight','bold','FontUnits','normalized',...
        'FontSize',gopt.fontsize/512);
end

hold off
%
% end function
function [movie_frame] = project_s2_partition(N,varargin)
%PROJECT_S2_PARTITION Use projection to illustrate an EQ partition of S^2
%
%Syntax
% [movie_frame] = project_s2_partition(N,options);
%
%Description
% PROJECT_S2_PARTITION(N) uses projection to illustrate the partition of
% the unit sphere S^2 into N regions.
%
% MOVIE_FRAME = PROJECT_S2_PARTITION(N) sets MOVIE_FRAME to be an array of
% movie frames for use with MOVIE. The movie frames will contain the region by
% region build-up of the illustration.
%
% PROJECT_S2_PARTITION(N,'offset','extra') uses experimental extra offsets.
% For more detail on partition options, see HELP PARTITION_OPTIONS.
%
% PROJECT_S2_PARTITION(N,options) also recognizes a number of illustration
% options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% PROJECT_S2_PARTITION(N,'fontsize',size)
% Font size used in titles (numeric, default 16).
%
% PROJECT_S2_PARTITION(N,'title','long')
% PROJECT_S2_PARTITION(N,'title','short')
% Use long or short titles (default 'long').
%
% PROJECT_S2_PARTITION(N,'proj','stereo')
% PROJECT_S2_PARTITION(N,'proj','eqarea')
% Use stereographic or equal area projection (default 'stereo').
%
% PROJECT_S2_PARTITION(N,'points','show')
% PROJECT_S2_PARTITION(N,'points','hide')
% Show or hide center points (default 'show').
%
% For more detail on illustration options, see HELP ILLUSTRATION_OPTIONS.
%
%Examples
%
% >> project_s2_partition(10)
% >> frames = project_s2_partition(9,'offset','extra','proj','eqarea')
%
% frames =
%
% 1x10 struct array with fields:
%     cdata
%     colormap
%
% >> project_s2_partition(99,'proj','eqarea','points','hide')
%
%See also
% MOVIE, PARTITION_OPTIONS, ILLUSTRATION_OPTIONS, SHOW_S2_PARTITION,
% PROJECT_S3_PARTITION

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

pdefault.extra_offset =  false;

popt = partition_options(pdefault, varargin{:});

gdefault.fontsize = 16;
gdefault.stereo =        true;
gdefault.show_title =    true;
gdefault.long_title =    true;
gdefault.show_points =   true;

gopt = illustration_options(gdefault, varargin{:});

dim = 2;

Phi = 2*pi*(0:1/40:1);
X = cos(Phi);
Y = sin(Phi);
if gopt.stereo
    r = 0;
else
    r = (area_of_sphere(dim)/volume_of_ball(dim)).^(1/dim);
end
plot(r*X,r*Y,'k')
axis equal;hold on
colormap jet
grid off
axis off

if gopt.show_title
    if gopt.long_title
        if gopt.stereo
            s = 'Stereographic';
        else
            s = 'Equal area';
        end
        if gopt.show_points
            pointstr = ', showing the center point of each region';
        else
            pointstr = '';
        end

        titlestr = sprintf(...
        '%s projection of recursive zonal equal area partition of {S^2}\ninto %d regions%s.',...
        s,N,pointstr);
    else
        titlestr = sprintf('EQ(2,%d)',N);
    end
    title(titlestr, ...
    'FontWeight','bold','FontUnits','normalized','FontSize',gopt.fontsize/512);
end

if nargout > 0
    movie_frame(1) = getframe(gcf);
end

R = eq_regions(dim,N,popt.extra_offset);
for i = N:-1:2
    project_s2_region(R(:,:,i),gopt.stereo);
    if nargout > 0
        movie_frame(N-i+2) = getframe(gcf);
    end
end

if gopt.show_points
    project_s2_eq_point_set(N,popt.extra_offset,gopt.stereo);
    if nargout > 0
        movie_frame(N+1) = getframe(gcf);
    end
end

hold off
%
% end function

function project_s2_region(region, stereo)
%PROJECT_S2_REGION Use projection to illustrate an EQ region of S^2
%
%Syntax
% project_s2_region(region, stereo);
%
%Notes
% The region is given as a 2 x 2 matrix in spherical polar coordinates
%
% The default is to use stereographic projection
% If the optional second argument, stereo is false,
% then use a equal area projection.

if nargin < 2
    stereo = true;
end
if stereo
    projection = @x2stereo;
else
    projection = @x2eqarea;
end

tol = eps*2^5;

dim = size(region,1);
t = region(:,1);
b = region(:,2);

if abs(b(1)) < tol
    b(1) = 2*pi;
end
pseudo = 0;
if abs(t(1)) < tol && abs(b(1)-2*pi) < tol
    pseudo = 1;
end
n = 21;
delta = 1/(n-1);
h = 0:delta:1;
t_to_b = zeros(dim,n);
b_to_t = t_to_b;
j = zeros(1,dim);
for k = 1:dim
    if ~pseudo || k < 2
        L = 1:dim;
        j(L) = mod(k+L,dim)+1;
        t_to_b(j(1),:) = t(j(1))+(b(j(1))-t(j(1)))*h;
        t_to_b(j(2),:) = t(j(2))*ones(1,n);
        t_to_b_x = polar2cart(t_to_b);
        s = projection(t_to_b_x);
        plot(s(1,:),s(2,:),'k');
        axis equal; hold on
        if pseudo
            axis equal; hold on
            b_to_t(j(1),:,:) = b(j(1))-(b(j(1))-t(j(1)))*h;
            b_to_t(j(2),:,:) = b(j(2))*ones(1,n);
            b_to_t_x = polar2cart(b_to_t);
            s = projection(b_to_t_x);
            plot(s(1,:),s(2,:),'k');
        end
    end
end
grid off
axis off
%
% end function

function project_s2_eq_point_set(N,extra_offset,stereo)
%PROJECT_S2_EQ_POINT_SET Use projection to illustrate an EQ point set of S^2
%
%Syntax
% project_s2_eq_point_set(N,extra_offset,stereo);

if nargin < 2
    extra_offset = true;
end
if nargin < 3
    stereo = true;
end
if stereo
    projection = 'stereo';
else
    projection = 'eqarea';
end

x = eq_point_set(2,N,extra_offset);
project_point_set(x,'title','hide','proj',projection);
hold on
%
% end function
function [movie_frame] = project_s3_partition(N,varargin)
%PROJECT_S3_PARTITION Use projection to illustrate an EQ partition of S^3
%
%Syntax
% [movie_frame] = project_s3_partition(N,options);
%
%Description
% PROJECT_S3_PARTITION(N) uses projection to illustrate the partition of
% the unit sphere S^3 into N regions.
%
% MOVIE_FRAME = PROJECT_S3_PARTITION(N) sets MOVIE_FRAME to be an array of
% movie frames for use with MOVIE. The movie frames will contain the region by
% region build-up of the illustration.
%
% PROJECT_S3_PARTITION(N,'offset','extra') uses experimental extra offsets.
% For more detail on partition options, see HELP PARTITION_OPTIONS.
%
% PROJECT_S3_PARTITION(N,options) also recognizes a number of illustration
% options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% PROJECT_S3_PARTITION(N,'fontsize',size)
% Font size used in titles (numeric, default 18).
%
% PROJECT_S3_PARTITION(N,'title','long')
% PROJECT_S3_PARTITION(N,'title','short')
% Use long or short titles (default 'long').
%
% PROJECT_S3_PARTITION(N,'proj','stereo')
% PROJECT_S3_PARTITION(N,'proj','eqarea')
% Use stereographic or equal area projection (default 'stereo').
%
% PROJECT_S3_PARTITION(N,'points','show')
% PROJECT_S3_PARTITION(N,'points','hide')
% Show or hide center points (default 'show').
%
% PROJECT_S3_PARTITION(N,'surf','show')
% PROJECT_S3_PARTITION(N,'surf','hide')
% Show or hide surfaces of regions (default 'show').
%
% For more detail on illustration options, see HELP ILLUSTRATION_OPTIONS.
%
%Examples
%
% >> project_s3_partition(10)
% >> frames = project_s3_partition(9,'offset','extra','proj','eqarea')
%
% frames =
%
% 1x18 struct array with fields:
%     cdata
%     colormap
%
% >> project_s3_partition(99,'proj','eqarea','points','hide')
%
%See also
% MOVIE, PARTITION_OPTIONS, ILLUSTRATION_OPTIONS, PROJECT_S2_PARTITION

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Function changed name from x2s2 to cart2polar2
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

pdefault.extra_offset =  false;

popt = partition_options(pdefault, varargin{:});

gdefault.fontsize = 18;
gdefault.stereo =        true;
gdefault.long_title =    true;
gdefault.show_points =   true;
gdefault.show_surfaces = true;

gopt = illustration_options(gdefault, varargin{:});

dim = 3;

[X,Y,Z] = sphere(90);
if gopt.stereo
    r = 0;
else
    r = (area_of_sphere(dim)/volume_of_ball(dim)).^(1/dim);
end

hold off

if gopt.show_surfaces
    surf(r*X,r*Y,r*Z,zeros(size(Z)),...
        'FaceAlpha',1/20,'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
else
    plot3(0,0,0,'w.')
end
axis equal;hold on
camlight right
colormap jet
grid off
axis off

if gopt.long_title
    if gopt.stereo
        s = 'Stereographic';
    else
        s = 'Equal volume';
    end
    if gopt.show_points
        pointstr = ', showing the center point of each region';
    else
        pointstr = '';
    end

    title(sprintf(...
        '\n%s projection of recursive zonal equal area partition of {S^3} \n into %d regions%s.',...
        s,N,pointstr),'FontSize',gopt.fontsize);
else
    title(sprintf('\nEQ(3,%d)',N),'FontSize',gopt.fontsize);
end

axis equal
grid off
axis off

pause(0);
if nargout > 0
    movie_frame(1) = getframe(gcf);
end

if gopt.stereo && (N == 1)
    return;
end

if popt.extra_offset
    [R,dim_1_rot] = eq_regions(dim,N,popt.extra_offset);
else
    R = eq_regions(dim,N);
end

for i = N:-1:2
    if popt.extra_offset
        project_s3_region(R(:,:,i),N,gopt.stereo,gopt.show_surfaces,dim_1_rot{i});
    else
        project_s3_region(R(:,:,i),N,gopt.stereo,gopt.show_surfaces);
    end
    pause(0);
    if nargout > 0
        movie_frame(N-i+2) = getframe(gcf);
    end
end

if gopt.show_points
    project_s3_eq_point_set(N,popt.extra_offset,gopt.stereo);
    if nargout > 0
        for k=1:min(N,40)
            movie_frame(N+k) = getframe(gcf);
        end
    end
end

hold off
%
% end function

function project_s3_region(region, N, stereo, show_surfaces, rot_matrix)
%PROJECT_S3_REGION Use projection to illustrate an EQ region of S^3
%Syntax
% project_s3_region(region, stereo, show_surfaces, rot_matrix);
%
%Notes
% The region is given as a 3 x 2 matrix in spherical polar coordinates
%
% The default is to use stereographic projection
% If the optional second argument, stereo is false,
% then use an equal area projection.

if nargin < 3
    stereo = true;
end
if stereo
    projection = @x2stereo;
else
    projection = @x2eqarea;
end
if nargin < 4
    show_surfaces = true;
end

offset_regions = (nargin >= 5);

tol = eps*2^5;

dim = size(region,1);
t = region(:,1);
b = region(:,2);

if abs(b(1)) < tol
    b(1) = 2*pi;
end
pseudo = 0;
if abs(t(1)) < tol && abs(b(1)-2*pi) < tol
    pseudo = 1;
end
n = 33;
delta = 1/(n-1);
h = 0:delta:1;
[h1, h2]  = meshgrid(h,h);
t_to_b = zeros(dim,n,n);
b_to_t = t_to_b;
r = N^(-1/3)/32;
L = 1:dim;
j = zeros(1,dim);
for k = 1:dim
    if ~pseudo || k < 3
        j(L) = mod(k+L,dim)+1;
        t_to_b(j(1),:,:) = t(j(1))+(b(j(1))-t(j(1)))*h1;
        t_to_b(j(2),:,:) = t(j(2))+(b(j(2))-t(j(2)))*h2;
        t_to_b(j(3),:,:) = t(j(3))*ones(n,n);
        t_to_b_v = reshape(t_to_b,dim,n*n);
        if offset_regions
            t_to_b_x = polar2cart([cart2polar2(rot_matrix*polar2cart(t_to_b_v(1:dim-1,:)));t_to_b_v(dim,:)]);
        else
            t_to_b_x = polar2cart(t_to_b_v);
        end
        s = reshape(projection(t_to_b_x),dim,n,n);
        degenerate = (norm(s(:,1,1)-s(:,1,2)) < tol);
        if ~degenerate && (~pseudo || k > 1)
            [X,Y,Z] = fatcurve(squeeze(s(:,1,:)),r);
            surface(X,Y,Z,zeros(size(Z)),...
                'FaceAlpha',(t(dim)/pi)/2,'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
            axis equal; hold on
        end
        if show_surfaces
            surf(squeeze(s(1,:,:)),squeeze(s(2,:,:)),squeeze(s(3,:,:)),t(3)*ones(n,n),...
                'FaceAlpha',(t(dim)/pi)/2,'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
        end
        axis equal; hold on

        b_to_t(j(1),:,:) = b(j(1))-(b(j(1))-t(j(1)))*h1;
        b_to_t(j(2),:,:) = b(j(2))-(b(j(2))-t(j(2)))*h2;
        b_to_t(j(3),:,:) = b(j(3))*ones(n,n);
        b_to_t_v = reshape(b_to_t,dim,n*n);
        if offset_regions
            b_to_t_x = polar2cart([cart2polar2(rot_matrix*polar2cart(b_to_t_v(1:dim-1,:)));b_to_t_v(dim,:)]);
        else
            b_to_t_x = polar2cart(b_to_t_v);
        end
        s = reshape(projection(b_to_t_x),dim,n,n);
        degenerate = (norm(s(:,1,1)-s(:,1,2)) < tol);
        if ~degenerate && (~pseudo || (k > 1 && abs(b(2)-pi) > tol))
            [X,Y,Z] = fatcurve(squeeze(s(:,1,:)),r);
            surface(X,Y,Z,zeros(size(Z)),...
                'FaceAlpha',(t(dim)/pi)/2,'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
        end
        if show_surfaces && k < 2
            surf(squeeze(s(1,:,:)),squeeze(s(2,:,:)),squeeze(s(3,:,:)),t(3)*ones(n,n),...
                'FaceAlpha',(t(dim)/pi)/2,'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
        end
    end
end
colormap jet
grid off
axis off
%
% end function

function project_s3_eq_point_set(N,extra_offset,stereo)
%PROJECT_S3_EQ_POINT_SET Use projection to illustrate an EQ point set of S^3
%
%Syntax
% project_s3_eq_point_set(N,min_energy,stereo);

if nargin < 2
    extra_offset = true;
end
if nargin < 3
    stereo = true;
end
if stereo
    projection = 'stereo';
else
    projection = 'eqarea';
end

x = eq_point_set(3,N,extra_offset);
project_point_set(x,'title','hide','proj',projection);
%
% end function
function show_r3_point_set(points_x,varargin)
%SHOW_R3_POINT_SET 3D illustration of a point set
%
%Syntax
% show_r3_point_set(POINTS_X,options);
%
%Description
% SHOW_R3_POINT_SET(POINTS_X) uses a 3d plot to illustrate a point set in relation to
% the unit sphere S^2.
%
% The argument POINTS_X must be an array of real numbers of size (3 by N), where N is a
% positive integer, representing N points of R^3.
%
% SHOW_R3_POINT_SET(POINTS_X,options) also recognizes a number of illustration
% options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% SHOW_R3_POINT_SET(POINTS_X,'fontsize',size)
% Font size used in titles (numeric, default 16).
%
% SHOW_R3_POINT_SET(POINTS_X,'title','show')
% SHOW_R3_POINT_SET(POINTS_X,'title','hide')
% Show or hide title (default 'hide').
%
% SHOW_R3_POINT_SET(POINTS_X,'sphere','show')
% SHOW_R3_POINT_SET(POINTS_X,'sphere','hide')
% Show or hide the unit sphere S^2 (default 'hide').
%
% For more detail on illustration options, see HELP ILLUSTRATION_OPTIONS.
%
%Note
% This function is primarily for use with the point set POINTS_X as a subset of the
% unit sphere S^2, but this is not assumed and not checked.
% If you show the unit sphere S^2 and POINTS_X contains points closer than radius 1
% from the origin, the sphere will hide these points.
%
%Examples
%
% >> points_x = [[0 0 1]' [0 1 0]' [0 -1 0]' [0 0 -1]']
%
% points_x =
%
%      0     0     0     0
%      0     1    -1     0
%      1     0     0    -1
%
% >> show_r3_point_set(points_x,'sphere','hide')
% >> show_r3_point_set(points_x,'sphere','show')
%
%See also
% ILLUSTRATION_OPTIONS, SHOW_S2_PARTITION, PROJECT_POINT_SET

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

gdefault.fontsize = 16;
gdefault.show_title  = false;
gdefault.show_sphere = false;

gopt = illustration_options(gdefault, varargin{:});

N = size(points_x,2);

surf_jet;

if gopt.show_title
    titlestr = sprintf(...
        '\nPoint set containing %d points.',N);
    title(titlestr,'FontWeight','bold','FontUnits','normalized',...
        'FontSize',gopt.fontsize/512);
end

if gopt.show_sphere
    show_s2_sphere;
    hold on
end

[X,Y,Z] = sphere;

r = min(0.05,N^(-1/2)/2);
rX = r*X; rY = r*Y; rZ = r*Z;
for n = 1:N
   surf(points_x(1,n)+rX,points_x(2,n)+rY,points_x(3,n)+rZ,ones(size(rZ)),...
   'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
end
%
axis equal
axis off
grid off
hold off
function [movie_frame] = show_s2_partition(N,varargin)
%SHOW_S2_PARTITION 3D illustration of an EQ partition of S^2
%
%Syntax
% [movie_frame] = show_s2_partition(N,options);
%
%Description
% SHOW_S2_PARTITION(N) uses a 3d plot to illustrate the partition of
% the unit sphere S^2 into N regions.
%
% MOVIE_FRAME = SHOW_S2_PARTITION(N) sets MOVIE_FRAME to be an array of
% movie frames for use with MOVIE. The movie frames will contain the region by
% region build-up of the illustration.
%
% SHOW_S2_PARTITION(N,'offset','extra') uses experimental extra offsets.
% For more detail on partition options, see HELP PARTITION_OPTIONS.
%
% SHOW_S2_PARTITION(N,options) also recognizes a number of illustration
% options, which are specified as name, value pairs.
% Any number of pairs can be used, in any order.
%
% The following illustration options are used.
%
% SHOW_S2_PARTITION(N,'fontsize',size)
% Font size used in titles (numeric, default 16).
%
% SHOW_S2_PARTITION(N,'title','show')
% SHOW_S2_PARTITION(N,'title','hide')
% Show or hide title (default 'show').
%
% SHOW_S2_PARTITION(N,'points','show')
% SHOW_S2_PARTITION(N,'points','hide')
% Show or hide center points (default 'show').
%
% SHOW_S2_PARTITION(N,'sphere','show')
% SHOW_S2_PARTITION(N,'sphere','hide')
% Show or hide the unit sphere S^2 (default 'show').
%
% For more detail on illustration options, see HELP ILLUSTRATION_OPTIONS.
%
%Examples
%
% >> show_s2_partition(10)
% >> frames = show_s2_partition(9,'offset','extra')
%
% frames =
%
% 1x10 struct array with fields:
%     cdata
%     colormap
%
% >> show_s2_partition(99,'points','hide')
%
%See also
% MOVIE, PARTITION_OPTIONS, ILLUSTRATION_OPTIONS, PROJECT_S2_PARTITION

% Copyright 2025 Paul Leopardi.
% $Revision 1.12.1 $ $Date 2025-08-09 $
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-16 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

pdefault.extra_offset =  false;
popt = partition_options(pdefault, varargin{:});

gdefault.fontsize = 16;
gdefault.show_title  = true;
gdefault.show_points = true;
gdefault.show_sphere = true;
gopt = illustration_options(gdefault, varargin{:});

dim = 2;

surf_jet;

if gopt.show_title
    if gopt.show_points
        pointstr = ', showing the center point of each region';
    else
        pointstr = '';
    end
    titlestr = sprintf(...
        '\nRecursive zonal equal area partition of {S^2} \n into %d regions%s.',...
        N,pointstr);
    title(titlestr,'FontWeight','bold','FontUnits','normalized',...
        'FontSize',gopt.fontsize/512);
end

frame_no = 1;
if gopt.show_sphere
    show_s2_sphere;
    hold on
    if nargout > 0
        movie_frame(frame_no) = getframe(gcf);
        frame_no = frame_no + 1;
    end
end

R = eq_regions(dim,N,popt.extra_offset);
top_colat = 0;
for i = N:-1:2
    if top_colat ~= R(2,1,i)
        top_colat = R(2,1,i);
        pause(0);
    end
    show_s2_region(R(:,:,i),N);
    if nargout > 0
        movie_frame(frame_no) = getframe(gcf);
        frame_no = frame_no + 1;
    end
end

if gopt.show_points
    x = eq_point_set(dim,N,popt.extra_offset);
    show_r3_point_set(x,'sphere','hide','title','hide');
    hold on
    if nargout > 0
        movie_frame(frame_no) = getframe(gcf);
    end
end

hold off
%
% end function

function show_s2_region(region,N)
%SHOW_S2_REGION Illustrate a region of S^2
%
%Syntax
% show_s2_region(region,N);
%
%Description
% SHOW_S2_REGION(REGION,N) uses 3D surface plots to illustrate a region of S^2.
% The region is given as a 2 x 2 matrix in spherical polar coordinates

tol = eps*2^5;

dim = size(region,1);
t = region(:,1);
b = region(:,2);

if abs(b(1)) < tol
    b(1) = 2*pi;
end
pseudo = 0;
if abs(t(1)) < tol && abs(b(1)-2*pi) < tol
    pseudo = 1;
end
n = 21;
delta = 1/(n-1);
h = 0:delta:1;
t_to_b = zeros(dim,n);
r = sqrt(1/N)/12;
j = zeros(1,dim);
for k = 1:dim
    if ~pseudo || k < 2
        L = 1:dim;
        j(L) = mod(k+L,dim)+1;
        t_to_b(j(1),:) = t(j(1))+(b(j(1))-t(j(1)))*h;
        t_to_b(j(2),:) = t(j(2))*ones(1,n);
        t_to_b_x = polar2cart(t_to_b);
        [X,Y,Z] = fatcurve(t_to_b_x,r);
        surface(X,Y,Z,-ones(size(Z)),...
       'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
        axis equal
        hold on
    end
end
grid off
axis off
%
% end function
