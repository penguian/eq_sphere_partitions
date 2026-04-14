function benchmark_sradius(varargin)
%BENCHMARK_SRADIUS Benchmark for sradius_of_cap performance.
%
% Usage:
%   benchmark_sradius('n_max', 1e7, 'dim', 3)

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
    addParameter(p, 'dim', 3);
    parse(p, varargin{:});

    args = p.Results;

    fprintf('%-10s | %10s\n', 'Size', 'Time (s)');
    fprintf('%s\n', repmat('-', 1, 23));

    n_values = generate_125_sequence(args.n_max);

    % Warm-up to absorb startup overhead
    for n_warm = [10, 20, 50]
        areas_warmup = linspace(0.1, area_of_sphere(args.dim) - 0.1, n_warm);
        sradius_of_cap(args.dim, areas_warmup);
    end

    results = [];
    for current_size = n_values
        areas = linspace(0.1, area_of_sphere(args.dim) - 0.1, current_size);

        t0 = tic;
        sradius_of_cap(args.dim, areas);
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
            fprintf('Best fitting order: O(Size^%.2f)\n', p_fit(1));
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
