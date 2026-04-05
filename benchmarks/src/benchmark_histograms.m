function benchmark_histograms(varargin)
%BENCHMARK_HISTOGRAMS Benchmark for point-to-region mapping (histogram).
%
% Usage:
%   benchmark_histograms('n_max', 2e8, 'regions', 1000)

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2026-04-05 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

    p = inputParser;
    addParameter(p, 'n_max', 2e8);
    addParameter(p, 'regions', 1000);
    parse(p, varargin{:});

    args = p.Results;

    fprintf('%-10s | %10s\n', 'Points', 'Time (s)');
    fprintf('%s\n', repmat('-', 1, 23));

    % Warm-up to absorb startup overhead
    test_points_warmup = [2*pi*rand(1, 100); pi*rand(1, 100)];
    eq_find_s2_region(test_points_warmup, args.regions);

    chunk_size = max(1, floor(args.n_max / 5));
    for current_size = max(100, chunk_size):chunk_size:args.n_max
        % Generate random points as test data on S^2
        test_points = [2*pi*rand(1, current_size); pi*rand(1, current_size)];

        t0 = tic;
        eq_find_s2_region(test_points, args.regions);
        t_elapsed = toc(t0);

        fprintf('%-10d | %10.4f\n', current_size, t_elapsed);
    end
end
