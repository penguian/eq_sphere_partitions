function benchmark_sradius(varargin)
%BENCHMARK_SRADIUS Benchmark for sradius_of_cap (Root finding loop bottleneck).
%
% Usage:
%   benchmark_sradius('n_max', 1e8, 'dim', 3)

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2026-04-05 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

    p = inputParser;
    addParameter(p, 'n_max', 1e8);
    addParameter(p, 'dim', 3);
    parse(p, varargin{:});

    args = p.Results;

    fprintf('%-10s | %10s\n', 'Size', 'Time (s)');
    fprintf('%s\n', repmat('-', 1, 23));

    % Warm-up to absorb startup overhead
    areas_warmup = linspace(0.1, area_of_sphere(args.dim) - 0.1, 100);
    sradius_of_cap(args.dim, areas_warmup);

    chunk_size = max(1, floor(args.n_max / 5));
    for current_size = max(100, chunk_size):chunk_size:args.n_max
        areas = linspace(0.1, area_of_sphere(args.dim) - 0.1, current_size);
        t0 = tic;
        sradius_of_cap(args.dim, areas);
        t_elapsed = toc(t0);

        fprintf('%-10d | %10.4f\n', current_size, t_elapsed);
    end
end
