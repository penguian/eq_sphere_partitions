function test_bot_cap_region
%TEST_BOT_CAP_REGION Test the bot_cap_region function

% Copyright 2025 Paul Leopardi.
% $Revision $ 1.12.2 $ $Date 2025-09-20 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

help bot_cap_region
region = bot_cap_region(1, pi/6)
region = bot_cap_region(2, pi/6)
region = bot_cap_region(3, pi/6)
