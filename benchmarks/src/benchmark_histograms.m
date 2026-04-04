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

    fprintf('Testing eq_find_s2_region: %d regions, %g points\n', ...
        args.regions, args.n_max);

    % Generate random points as test data on S^2
    % eq_find_s2_region expects (2 x n_points) in polar coordinates
    n_test = 1e6;
    fprintf('Generating %d test points...\n', n_test);
    test_points = [2*pi*rand(1, n_test); pi*rand(1, n_test)];

    % Warm-up
    eq_find_s2_region([0; 0], args.regions);

    t0 = tic;
    r_idx = eq_find_s2_region(test_points, args.regions);
    t_elapsed = toc(t0);

    % Scaling to n_max
    total_time = t_elapsed * (args.n_max / n_test);

    fprintf('Processed %d points in %.4f seconds (scaled to %g points: %.4f seconds)\n', ...
        n_test, t_elapsed, args.n_max, total_time);
    fprintf('Throughput: %g points/sec\n', n_test / t_elapsed);

    % Verify result size
    if length(r_idx) ~= n_test
        warning('Unexpected result size: %d', length(r_idx));
    end
end
