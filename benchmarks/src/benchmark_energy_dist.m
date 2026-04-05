function benchmark_energy_dist(varargin)
%BENCHMARK_ENERGY_DIST Benchmark for point_set_energy_dist (eq_energy_dist).
%
% Usage:
%   benchmark_energy_dist('n_max', 2400, 'dim', 2, 's', 2, 'extra_offset', false)

    p = inputParser;
    addParameter(p, 'n_max', 2400);
    addParameter(p, 'dim', 2);
    addParameter(p, 's', 2);
    addParameter(p, 'extra_offset', false);
    parse(p, varargin{:});
    
    args = p.Results;
    
    fprintf('%-15s | %10s\n', 'N range', 'Time (s)');
    fprintf('%s\n', repmat('-', 1, 28));
    
    % Warm-up from N=10 to 99 to absorb startup overhead
    for N = 10:99
        eq_energy_dist(args.dim, N, args.s, args.extra_offset);
    end

    chunk_size = max(1, floor(args.n_max / 5));
    
    for i = 0:chunk_size:args.n_max-1
        start_n = max(100, i + 1);
        end_n = min(i + chunk_size, args.n_max);
        
        if start_n > end_n
            continue;
        end
        
        t0 = tic;
        for N = start_n:end_n
            eq_energy_dist(args.dim, N, args.s, args.extra_offset);
        end
        t_elapsed = toc(t0);
        
        fprintf('%-15s | %10.4f\n', sprintf('%d-%d', start_n, end_n), t_elapsed);
    end
end
