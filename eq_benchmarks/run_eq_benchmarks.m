function run_benchmarks(varargin)
%RUN_BENCHMARKS Run all performance benchmarks and log results.
%
% This script runs a subset of benchmark points using logarithmic sampling
% (powers of 10) to verify asymptotic behavior and find the runtime order.
%
% Usage:
%   run_benchmarks('n_max_regions', 1e8)

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.3 $ $Date 2026-04-05 $
%
% For licensing, see COPYING.

    p = inputParser;
    addParameter(p, 'dim', 2);
    addParameter(p, 'n_max_regions', 1e8);
    addParameter(p, 'n_max_hist', 2e6);
    addParameter(p, 'n_max_srad', 1e6);
    addParameter(p, 'n_max_area', 5e7);
    addParameter(p, 'n_max_energy', 5e4);
    addParameter(p, 'n_max_mindist', 5e4);
    addParameter(p, 'extra_offset', false);
    parse(p, varargin{:});
    args = p.Results;

    % Ensure results directory exists
    if ~exist('benchmarks/results', 'dir')
        mkdir('benchmarks/results');
    end

    main_log_file = fullfile('benchmarks', 'results', 'run_benchmarks.log');

    % Start logging to file (diary)
    if exist(main_log_file, 'file')
        delete(main_log_file);
    end
    diary(main_log_file);

    fprintf('=======================================\n');
    fprintf('    MATLAB EQSP Performance Benchmarks \n');
    fprintf('=======================================\n\n');

    fprintf('Hardware: AMD Ryzen 7 8840HS w/ Radeon 780M Graphics (~2.4 GHz)\n');
    fprintf('OS:       Linux\n');
    fprintf('Software: MATLAB R2023b\n\n');

    t_overall_start = tic;

    % 1. eq_area_error
    fprintf('\nRunning benchmark: eq_area_error\n');
    n_values = generate_sequence(args.n_max_area);
    run_sampled_task(@(N) benchmark_area_error_point(args.dim, N, args.extra_offset), n_values, 'N');

    % 2. eq_regions
    fprintf('\nRunning benchmark: eq_regions\n');
    n_values = generate_sequence(args.n_max_regions);
    run_sampled_task(@(N) eq_regions(args.dim, N, args.extra_offset), n_values, 'N');

    % 3. eq_find_s2_region
    fprintf('\nRunning benchmark: eq_find_s2_region\n');
    n_values = generate_sequence(args.n_max_hist);
    run_sampled_task(@(N) benchmark_histograms_point(N, 1000), n_values, 'Points');

    % 4. sradius_of_cap
    fprintf('\nRunning benchmark: sradius_of_cap\n');
    n_values = generate_sequence(args.n_max_srad);
    run_sampled_task(@(N) benchmark_sradius_point(args.dim + 1, N), n_values, 'Size');

    % 5. eq_min_dist
    fprintf('\nRunning benchmark: eq_min_dist\n');
    n_values = generate_sequence(args.n_max_mindist);
    run_sampled_task(@(N) eq_min_dist(args.dim, N, args.extra_offset), n_values, 'N');

    % 6. point_set_energy_dist
    fprintf('\nRunning benchmark: point_set_energy_dist\n');
    n_values = generate_sequence(args.n_max_energy);
    run_sampled_task(@(N) eq_energy_dist(args.dim, N, 2, args.extra_offset), n_values, 'N');

    t_overall_end = toc(t_overall_start);
    fprintf('\n=======================================\n');
    fprintf('Total benchmark time: %.2f seconds\n', t_overall_end);
    fprintf('=======================================\n');

    diary off;
    fprintf('\nResults saved to %s\n', main_log_file);
end

function run_sampled_task(task_func, n_values, label)
    fprintf('%-15s | %10s\n', label, 'Time (s)');
    fprintf('%s\n', repmat('-', 1, 28));
    
    results = zeros(length(n_values), 2);
    for i = 1:length(n_values)
        N = n_values(i);
        t0 = tic;
        task_func(N);
        t_elapsed = toc(t0);
        results(i, :) = [N, t_elapsed];
        if N >= 100
            fprintf('%-15d | %10.4f\n', N, t_elapsed);
        end
    end
    
    % Scaling Analysis: Find best fitting power x in O(N^x)
    % Use only points where N >= 100 to avoid cold-start overhead, if available
    mask = (results(:,1) >= 100) & (results(:,2) > 0.0001);
    if sum(mask) < 2
        mask = (results(:,2) > 0); % Fallback to all non-zero points
    end
    
    if sum(mask) >= 2
        log_n = log(results(mask, 1));
        log_t = log(results(mask, 2));
        p = polyfit(log_n, log_t, 1);
        x = p(1);
        fprintf('Best fitting order: O(N^%.2f)\n', x);
    else
        fprintf('Insufficient data for scaling analysis.\n');
    end
end

function n_values = generate_sequence(n_max)
    % Generate 1-2-5 sequence: 5, 10, 20, 50, 100, 200, 500, ...
    v = [1, 2, 5]';
    p = 0:ceil(log10(n_max));
    vals = v * 10.^p;
    n_values = sort(vals(:))';
    n_values = n_values(n_values >= 10 & n_values <= n_max);
end

function benchmark_area_error_point(dim, N, extra_offset)
    regions = eq_regions(dim, N, extra_offset);
    % Simulate some work with the regions to avoid optimization skipping
    size(regions);
end

function benchmark_histograms_point(n_points, regions)
    test_points = [2*pi*rand(1, n_points); pi*rand(1, n_points)];
    eq_find_s2_region(test_points, regions);
end

function benchmark_sradius_point(dim, n_areas)
    areas = linspace(0.1, area_of_sphere(dim) - 0.1, n_areas);
    sradius_of_cap(dim, areas);
end
