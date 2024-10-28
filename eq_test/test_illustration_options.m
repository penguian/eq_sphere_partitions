function test_illustration_options
%TEST_ILLUSTRATION_OPTIONS Test the illustration_options function

% Copyright 2024 Paul Leopardi.
% $Revision $ 1.12 $ $Date 2024-10-14 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help illustration_options
gdefault.fontsize=14;
gopt = illustration_options(gdefault,'proj','stereo')
gopt = illustration_options(gdefault,'proj','stereo','fontsize',12)
