function test_partition_options
%TEST_PARTITION_OPTIONS Test the partition_options function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help partition_options
pdefault.extra_offset = false;
popt = partition_options(pdefault,'offset','extra')
popt = partition_options(pdefault,false)
