function benchmark_eq_regions(varargin)
%BENCHMARK_EQ_REGIONS Replicate Thesis Benchmark for eq_regions.
%
% This script measures the performance of the eq_regions function,
% following the parameters in the pyeqsp benchmark.
%
% Usage:
%   benchmark_eq_regions('max_d', 4, 'max_k', 18, 'iterations', 3, 'show_progress', true)

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2026-04-05 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

    p = inputParser;
    addParameter(p, 'max_d', 4);
    addParameter(p, 'max_k', 18);
    addParameter(p, 'iterations', 3);
    addParameter(p, 'show_progress', false);
    addParameter(p, 'extra_offset', false);
    parse(p, varargin{:});

    args = p.Results;

    fprintf('Running Thesis Benchmark for eq_regions: d=[1..%d], N=2^[1..%d]\n', ...
        args.max_d, args.max_k);
    fprintf('Iterations per point: %d\n', args.iterations);

    if args.show_progress
        fprintf('%-3s | %-3s | %-10s | %15s\n', 'd', 'k', 'N', 'Mean Time (s)');
        fprintf('%s\n', repmat('-', 1, 45));
    end

    % Warm-up to let JIT compile the functions
    if args.show_progress
        fprintf('Warming up...\n');
    end
    for d = 1:args.max_d
        eq_regions(d, 2^8, args.extra_offset);
    end

    results = [];

    start_total = tic;

    for d = 1:args.max_d
        for k = 1:args.max_k
            N = 2^k;
            times = zeros(1, args.iterations);
            for i = 1:args.iterations
                t0 = tic;
                eq_regions(d, N, args.extra_offset);
                times(i) = toc(t0);
            end

            mean_time = mean(times);
            results = [results; d, k, N, mean_time]; %#ok<AGROW>

            if args.show_progress
                fprintf('%-3d | %-3d | %-10d | %15.6f\n', d, k, N, mean_time);
            end
        end
    end

    end_total = toc(start_total);

    % Analyze scaling for k >= 10 to capture asymptotic behavior
    analyze_scaling(results(results(:,2) >= 10, :));

    fprintf('\nTotal benchmark wall time: %.2f seconds\n', end_total);
end

function analyze_scaling(results)
    fprintf('\nScaling Analysis (t ~ N^x) for eq_regions:\n');

    dims = unique(results(:,1));
    for i = 1:length(dims)
        d = dims(i);
        d_mask = (results(:,1) == d);
        d_results = results(d_mask, :);

        if size(d_results, 1) < 2
            continue;
        end

        log_n = log(d_results(:,3));
        log_t = log(d_results(:,4));

        % Linear regression for x: log(t) = log(C) + x * log(N)
        % polyfit(X, Y, 1) returns [slope, intercept]
        p = polyfit(log_n, log_t, 1);
        x = p(1);

        fprintf('Dimension %-2d: x = %.4f (Thesis baseline: ~0.60)\n', d, x);
    end
end
