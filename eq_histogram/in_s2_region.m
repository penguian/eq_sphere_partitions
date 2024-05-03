function result=in_s2_region(s_point,region)

% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-04-30 $
% Copyright 2012 Paul Leopardi
% $Revision 1.11 $ $Date 2012-01-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

n_points=size(s_point,2);
result=zeros(1,n_points);
for p_idx=1:n_points
 result(p_idx)=false;
 longitude=s_point(1,p_idx);
 min_long=region(1,1);
 max_long=region(1,2);
 if min_long < longitude && longitude <= max_long
  result(p_idx)=true;
 else
  longitude=longitude + 2*pi;
  if min_long < longitude && longitude <= max_long
   result(p_idx)=true;
  end
 end
 colatitude=s_point(2,p_idx);
 min_colat=region(2,1);
 max_colat=region(2,2);
 if min_colat >= colatitude || colatitude > max_colat
   result(p_idx)=false;
 end
end


