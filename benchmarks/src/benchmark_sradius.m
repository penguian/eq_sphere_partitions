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

    fprintf('Testing sradius_of_cap: dim=%d, %g calls\n', args.dim, args.n_max);

    n_test = 10000;

    % Warm-up
    sradius_of_cap(args.dim, 1/100);

    t0 = tic;
    for i = 1:n_test
        sradius_of_cap(args.dim, 1/i);
    end
    t_elapsed = toc(t0);

    % Scaling to n_max
    total_time = t_elapsed * (args.n_max / n_test);

    fprintf('Processed %d calls in %.4f seconds (scaled to %g calls: %.4f seconds)\n', ...
        n_test, t_elapsed, args.n_max, total_time);
    fprintf('Throughput: %g calls/sec\n', n_test / t_elapsed);
end
