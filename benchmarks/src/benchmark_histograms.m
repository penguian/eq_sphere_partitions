function benchmark_histograms(varargin)
%BENCHMARK_HISTOGRAMS Benchmark for point-to-region mapping (histogram).
%
% Usage:
%   benchmark_histograms('n_max', 1e7, 'regions', 1000)

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.3 $ $Date 2026-04-14 $
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

    p = inputParser;
    addParameter(p, 'n_max', 1e7);
    addParameter(p, 'regions', 1000);
    parse(p, varargin{:});

    args = p.Results;

    fprintf('%-10s | %10s\n', 'Points', 'Time (s)');
    fprintf('%s\n', repmat('-', 1, 23));

    n_values = generate_125_sequence(args.n_max);

    % Warm-up to absorb startup overhead
    for n_warm = [10, 20, 50]
        test_points_warmup = [2*pi*rand(1, n_warm); pi*rand(1, n_warm)];
        eq_find_s2_region(test_points_warmup, args.regions);
    end

    results = [];
    for current_size = n_values
        % Generate random points as test data on S^2
        test_points = [2*pi*rand(1, current_size); pi*rand(1, current_size)];

        t0 = tic;
        eq_find_s2_region(test_points, args.regions);
        t_elapsed = toc(t0);

        results = [results; current_size, t_elapsed];
        if current_size >= 100
            fprintf('%-10d | %10.4f\n', current_size, t_elapsed);
        end
    end

    % Scaling Analysis
    if size(results, 1) >= 2
        mask = (results(:, 1) >= 100) & (results(:, 2) > 0.0001);
        if sum(mask) >= 2
            log_n = log(results(mask, 1));
            log_t = log(results(mask, 2));
            p_fit = polyfit(log_n, log_t, 1);
            fprintf('%s\n', repmat('-', 1, 23));
            fprintf('Best fitting order: O(Points^%.2f)\n', p_fit(1));
        end
    end
end

function vals = generate_125_sequence(n_max)
    % Generate 1-2-5 logarithmic sequence: 10, 20, 50, 100...
    if n_max < 10
        vals = max(1, floor(n_max));
        return;
    end
    v = [1, 2, 5];
    p = 0:ceil(log10(n_max));
    [V, P] = meshgrid(v, p);
    vals = sort(reshape(V .* (10.^P), 1, []));
    vals = vals(vals >= 10 & vals <= n_max);
end
