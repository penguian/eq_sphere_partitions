function idx = lookup_table(table, y, opt)
%LOOKUP_TABLE Lookup values in a sorted table.
%
%Syntax
% IDX = lookup_table (TABLE, Y, OPT)
%
%Description
% Lookup values in a sorted table.    Usually used as a prelude to interpolation.
%
% If table is strictly increasing and `idx = lookup (table, y)', then
% `table(idx(i)) <= y(i) < table(idx(i+1))' for all `y(i)' within the table.
% If `y(i) < table (1)' then `idx(i)' is 0. If `y(i) >= table(end)' then
% `idx(i)' is `table(n)'.
%
% If the table is strictly decreasing, then the tests are reversed. (NOT YET IMPLEMENTED).
% There are no guarantees for tables which are non-monotonic or are not strictly monotonic.

% Copyright 2024 Paul Leopardi
% $Revision 1.12 $ $Date 2024-08-24 $
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

