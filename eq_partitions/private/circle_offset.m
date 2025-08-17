function offset = circle_offset(n_top,n_bot,extra_twist)
%CIRCLE_OFFSET Try to maximize minimum distance of center points for S^2 collars
%
%Syntax
% offset = circle_offset(n_top,n_bot,extra_twist);
%
%Description
% OFFSET = CIRCLE_OFFSET(N_TOP,N_BOT,EXTRA_TWIST) calculates the offset OFFSET
% that maximizes the minimum distance between two sets of points on a circle,
% determined by the numbers N_TOP and N_BOT.
%
%Notes
% The values N_TOP and N_BOT represent the numbers of points in each of two
% equally speced sets on the unit circle. The offset is given in multiples of
% whole rotations, and consists of three parts:
% 1) Half the difference between a twist of one sector on each of bottom and top.
% This brings the centre points into alignment.
% 2) A rotation which will maximize the minimum angle between points on the two
% circles.
% 3) An optional extra twist by a whole number of sectors on the second circle.
% The extra twist is added so that the location of the minimum angle between
% circles will progressively twist around the sphere with each collar.
%
%Examples
%
% >> offset = circle_offset(3, 4)
%
% offset =
%
%    6.9389e-18
%
% >> offset = circle_offset(3, 4, 1)
%
% offset =
%
%     1.5000
%
% >> offset = circle_offset(7, 11)
%
% offset =
%
%    -0.0195

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

if nargin < 3
    extra_twist = false;
end
offset = (1/n_bot - 1/n_top)/2 + gcd(n_top,n_bot)/(2*n_top*n_bot);
if extra_twist
    twist = 6;
    offset = offset + twist/n_bot;
end
%
% end function
