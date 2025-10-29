function n_collars = num_collars(N,c_polar,a_ideal)
%NUM_COLLARS The number of collars between the polar caps
%
%Syntax
%  n_collars = num_collars(N,c_polar,a_ideal);
%
%Description
% N_COLLARS = NUM_COLLARS(N,C_POLAR,A_IDEAL) determines N_COLLARS, the number
% of collars between the polar caps of an equal area partition of the sphere
% into N regions, given the North polar colatitude C_POLAR, and the ideal
% collar angle A_IDEAL.
%
%Examples
%
% >> dim = 2; N = 4;
% >> c_polar = polar_colat(dim, N);
% >> a_ideal = ideal_collar_angle(dim, N);
% >> n_collars = num_collars(N,c_polar,a_ideal)
%
% n_collars =
%
%      1
%
% >> dim = 2; N = 10;
% >> c_polar = polar_colat(dim, N);
% >> a_ideal = ideal_collar_angle(dim, N);
% >> n_collars = num_collars(N,c_polar,a_ideal)
%
% n_collars =
%
%      2
%
% >> dim = 3; N = 6;
% >> c_polar = polar_colat(dim, N);
% >> a_ideal = ideal_collar_angle(dim, N);
% >> n_collars = num_collars(N,c_polar,a_ideal)
%
% n_collars =
%
%      1
%
%See also
% IDEAL_COLLAR_ANGLE, IDEAL_REGION_LIST, POLAR_COLAT

% Copyright 2025 Paul Leopardi
% $Revision 1.12.2 $ $Date 2025-08-17 $
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_collars = zeros(size(N));
enough = (N > 2) & (a_ideal > 0);
n_collars(enough) = max(1,round((pi-2*c_polar(enough))./a_ideal(enough)));
%
% end function
