function benchmark_eq_regions(varargin)
%BENCHMARK_EQ_REGIONS Performance benchmark for eq_regions.
%
% This script measures the performance of the eq_regions function,
% following the parameters in the pyeqsp benchmark.
%
% Usage:
%   benchmark_eq_regions('n_max', 1e6, 'dim', 2)

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.3 $ $Date 2026-04-14 $
% Copyright 2024 Paul Leopardi.
% $Revision 1.12 $ $Date 2024-10-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

    p = inputParser;
    addParameter(p, 'n_max', 1e6);
    addParameter(p, 'dim', 2);
    parse(p, varargin{:});

    args = p.Results;

    fprintf('%-10s | %10s\n', 'N', 'Time (s)');
    fprintf('%s\n', repmat('-', 1, 23));

    n_values = generate_125_sequence(args.n_max);

    % Warm-up to absorb startup overhead
    for n_warm = [10, 20, 50]
        eq_regions(args.dim, n_warm);
    end

    results = zeros(numel(n_values), 2);
    result_idx = 0;
    for N = n_values
        t0 = tic;
        eq_regions(args.dim, N);
        t_elapsed = toc(t0);

        result_idx = result_idx + 1;
        results(result_idx, :) = [N, t_elapsed];
        if N >= 100
            fprintf('%-10d | %10.4f\n', N, t_elapsed);
        end
    end

    % Scaling Analysis
    if size(results, 1) >= 2
        mask = (results(:,1) >= 100) & (results(:,2) > 0.0001);
        if sum(mask) >= 2
            log_n = log(results(mask, 1));
            log_t = log(results(mask, 2));
            p_fit = polyfit(log_n, log_t, 1);
            fprintf('%s\n', repmat('-', 1, 23));
            fprintf('Best fitting order: O(N^%.2f)\n', p_fit(1));
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
