function r_idx = eq_find_s2_region(s_point,N)

% Copyright 2012 Paul Leopardi
% $Revision 1.11 $ $Date 2012-01-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

[s_cap, n_regions] = eq_caps(2,N);
s_regions = eq_regions(2,N);
n_caps = length(s_cap);
c_regions = n_regions;
for c_idx = 2:n_caps
 c_regions(c_idx) = c_regions(c_idx-1) + n_regions(c_idx);
end
r_idx = lookup_s2_region(s_point,s_regions,s_cap,c_regions);
