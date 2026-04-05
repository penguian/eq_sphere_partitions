function benchmark_mindist(varargin)
%BENCHMARK_MINDIST Benchmark for eq_min_dist (Efficient distance calculation).
%
% Usage:
%   benchmark_mindist('n_max', 16000, 'dim', 2, 'stepsize', 1, 'extra_offset', false)

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2026-04-05 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

    p = inputParser;
    addParameter(p, 'n_max', 6400);
    addParameter(p, 'dim', 2);
    addParameter(p, 'stepsize', 1);
    addParameter(p, 'extra_offset', false);
    parse(p, varargin{:});

    args = p.Results;

    fprintf('Testing eq_min_dist: dim=%d, N=1..%d (step %d)\n', args.dim, args.n_max, args.stepsize);
    fprintf('%-15s | %10s\n', 'N range', 'Time (s)');
    fprintf('%s\n', repmat('-', 1, 28));

    % Warm-up from N=10 to 99 to absorb startup overhead
    for N = 10:99
        eq_min_dist(args.dim, N, args.extra_offset);
    end

    chunk_size = max(1, floor(args.n_max / 5));

    for i = 0:chunk_size:args.n_max-1
        start_n = max(100, i + 1);
        end_n = min(i + chunk_size, args.n_max);

        if start_n > end_n
            continue;
        end

        t0 = tic;
        for N = start_n:args.stepsize:end_n
            eq_min_dist(args.dim, N, args.extra_offset);
        end
        t_elapsed = toc(t0);

        fprintf('%-15s | %10.4f\n', sprintf('%d-%d', start_n, end_n), t_elapsed);
    end
end
