function test_centres_of_regions
%TEST_CENTRES_OF_REGIONS Test the centres_of_regions function

% Copyright 2025 Paul Leopardi.
% $Revision $ 1.12.2 $ $Date 2025-08-17 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help centres_of_regions
dim = 2; N = 4;
regions = eq_regions(dim, N);
points = centres_of_regions(regions)
dim = 2; N = 10;
regions = eq_regions(dim, N);
points = centres_of_regions(regions)
dim = 3; N = 6;
regions = eq_regions(dim, N);
points = centres_of_regions(regions)
