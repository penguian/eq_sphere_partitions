function idx = lookup_table(table, y, opt)
%LOOKUP_TABLE Lookup values in a sorted table.
%
%Syntax
% idx = lookup_table(table, y, opt)
%
%Description
% Lookup values in a sorted table. Usually used as a prelude to interpolation.
%
% If table is strictly increasing and `idx = lookup (table, y)', then
% `table(idx(i)) <= y(i) < table(idx(i)+1)' for all `y(i)' within the interval 
% with minimum table(1) and maximum table(end). If `y(i) < table(1)' then 
% `idx(i)' is 0. If `y(i) >= table(end)' then `idx(i)' is `end'.
%
% If the table is strictly decreasing, then the tests are reversed. (NOT YET IMPLEMENTED).
% There are no guarantees for tables which are non-monotonic or are not strictly monotonic.
% 
%Arguments
% table    Sequence of real values, assumed to be strictly increasing or decreasing.
% y        Sequence of real values to be looked up in table.
% opt      Options to be passed to Octave lookup, if we are in Octave, otherwise ignored.
%
%Examples
%  > table = [-100.0, -70, 2.5, 75, 125.7]
%  table =
%
%    -100.0000  -70.0000    2.5000   75.0000  125.7000
%
%  > y = [-1, 3, 1000, -197]
%  y =
%
%           -1           3        1000        -197
%
%  > idx = lookup_table(table, y)
%  idx =
%
%       2     3     5     0

% Copyright 2024 Paul Leopardi
% $Revision 1.12 $ $Date 2024-09-06 $
% Copyright 2013 Paul Leopardi
% $Revision 1.11 $ $Date 2013-04-18 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

if is_octave
    if nargin == 3
        idx = lookup(table, y, opt);
    else
        idx = lookup(table, y);
    end
else
    if table(end) >= table(1)
        % We have an nondecreasing table.
        maximum = max([table, y]) + 1;
        idx = arrayfun(@(x) find(x < [table, maximum], 1) - 1, y);
    else
        fprintf('LOOKUP_TABLE: Decreasing case is not yet implemented\n');
    end
end

% Subfunction that checks if we are in Octave
function r = is_octave ()
    persistent x;
    if (isempty (x))
        x = exist ('OCTAVE_VERSION', 'builtin');
    end
    r = x;
    end
end

