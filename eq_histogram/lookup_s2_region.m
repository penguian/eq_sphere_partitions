function r_idx=lookup_s2_region(s_point,s_regions,s_cap,c_regions)

% Copyright 2012 Paul Leopardi
% $Revision 1.11 $ $Date 2012-01-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_caps=length(s_cap);
if n_caps ~= length(c_regions)
 fprintf('LOOKUP_S2_REGION: Mismatch between length of s_cap (==%d) and length of c_regions (==%d)\n',n_caps,length(c_regions))
 r_idx=0;
 return
end
n_regions=size(s_regions,3);
if c_regions(n_caps) ~= n_regions
 fprintf('LOOKUP_S2_REGION: Mismatch between c_regions(end) (==%d) and length of s_regions (==%d)\n',c_regions(n_caps),n_regions)
 r_idx=0;
 return
end
n_points=size(s_point,2);
r_idx=zeros(1,n_points);
for p_idx=1:n_points
 c_idx=mylookup(s_cap,s_point(2,p_idx));
 if c_idx>0 & c_idx<n_caps-1
  min_r_idx=c_regions(c_idx)+1;
  max_r_idx=c_regions(c_idx+1);
  s_longs=squeeze(s_regions(1,:,min_r_idx:max_r_idx));
  if s_longs(2,1) >= 2*pi
    s_longs(:,1) = s_longs(:,1) - 2*pi;
  end
  n_longs=size(s_longs,2);
  l_idx=mod(mylookup(s_longs(2,:),s_point(1,p_idx)),n_longs);
  if s_point(1,p_idx) < s_longs(1,1)
   l_idx=n_longs-1;
  end
  r_idx(p_idx)=min_r_idx+l_idx;
 elseif c_idx==0
  r_idx(p_idx)=1;
 elseif c_idx==n_caps-1
  r_idx(p_idx)=n_regions;
 else
  r_idx(p_idx)=0;
 end
end

