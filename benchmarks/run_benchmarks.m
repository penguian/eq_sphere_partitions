function run_benchmarks(varargin)
%RUN_BENCHMARKS Run all performance benchmarks and log results.
%
% Usage:
%   run_benchmarks('dim', 2, 'extra_offset', false)

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.2 $ $Date 2026-04-05 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

    p = inputParser;
    addParameter(p, 'dim', 2);
    addParameter(p, 'extra_offset', false);
    addParameter(p, 'mindist_n_max', 16000);
    addParameter(p, 'mindist_stepsize', 1);
    parse(p, varargin{:});
    args = p.Results;

    % Ensure directories exist
    if ~exist('benchmarks/results', 'dir')
        mkdir('benchmarks/results');
    end

    suffix = '';
    if args.extra_offset
        suffix = '_extra';
    end

    main_log_file = fullfile('benchmarks', 'results', ['run_benchmarks', suffix, '.log']);

    % Start logging to file (diary)
    if exist(main_log_file, 'file')
        delete(main_log_file);
    end
    diary(main_log_file);

    fprintf('=======================================\n');
    fprintf('     EQSP Performance Benchmarks       \n');
    fprintf('=======================================\n\n');

    t_start = tic;

    fprintf('\n--- Running eq_regions benchmark ---\n');
    benchmark_eq_regions('max_d', args.dim, 'max_k', 18, 'show_progress', true, 'extra_offset', args.extra_offset);

    fprintf('\n--- Running eq_area_error benchmark ---\n');
    benchmark_area_error('n_max', 15000, 'dim', args.dim, 'extra_offset', args.extra_offset);

    fprintf('\n--- Running eq_energy_dist benchmark ---\n');
    benchmark_energy_dist('n_max', 2400, 'dim', args.dim, 'extra_offset', args.extra_offset);

    fprintf('\n--- Running eq_min_dist benchmark ---\n');
    benchmark_mindist('n_max', args.mindist_n_max, 'dim', args.dim, 'stepsize', args.mindist_stepsize, 'extra_offset', args.extra_offset);

    fprintf('\n--- Running eq_histograms benchmark ---\n');
    benchmark_histograms('n_max', 2e8, 'regions', 1000);

    fprintf('\n--- Running sradius benchmark ---\n');
    benchmark_sradius('n_max', 1e8, 'dim', args.dim + 1);

    t_end = toc(t_start);
    fprintf('\n=======================================\n');
    fprintf('Total benchmark time: %.2f seconds\n', t_end);
    fprintf('=======================================\n');

    diary off;
    fprintf('\nResults saved to %s\n', main_log_file);
end
