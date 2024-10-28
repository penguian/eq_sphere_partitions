function test_project_s3_partition
%TEST_PROJECT_S3_PARTITION Test the project_s3_partition function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help project_s3_partition
project_s3_partition(10)
frames = project_s3_partition(9,'offset','extra','proj','eqarea')
project_s3_partition(99,'proj','eqarea','points','hide')
