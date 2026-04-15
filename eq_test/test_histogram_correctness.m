function test_histogram_correctness
%TEST_HISTOGRAM_CORRECTNESS Comprehensive tests for histogram functions
%
% Ported from pyeqsp tests:
% - test_consistency_find_and_in_region
% - test_boundary_conditions
% - test_exact_boundaries_s2_region

% Copyright 2026 Paul Leopardi.
% $Revision 1.12.3 $ $Date 2026-04-14 $

test_consistency;
test_poles;
test_boundaries;

fprintf('test_histogram_correctness: PASS\n');
end

function test_consistency
    N = 10;
    num_points = 50;
    % Generate random points on S^2 (seeded for determinism)
    rng_state = rng;
    rng(42);
    points_s = zeros(2, num_points);
    points_s(1, :) = rand(1, num_points) * 2 * pi; % Longitude
    points_s(2, :) = acos(rand(1, num_points) * 2 - 1); % Colatitude
    rng(rng_state);

    region_indices = eq_find_s2_region(points_s, N);
    regions = eq_regions(2, N);

    for i = 1:num_points
        r_idx = region_indices(i);
        region = regions(:, :, r_idx);
        assert(in_s2_region(points_s(:, i), region), ...
            sprintf('Point %d assigned to region %d but in_s2_region mismatch', i, r_idx));

        % Check a distant region
        other_idx = mod(r_idx + floor(N/2) - 1, N) + 1;
        other_region = regions(:, :, other_idx);
        assert(~in_s2_region(points_s(:, i), other_region), ...
            sprintf('Point %d assigned to region %d but also claims to be in region %d', i, r_idx, other_idx));
    end
end

function test_poles
    N = 4;
    % North Pole
    points_np = [0; 0];
    r_idx = eq_find_s2_region(points_np, N);
    assert(r_idx == 1, 'North pole should be in region 1');

    % South Pole
    points_sp = [0; pi];
    r_idx = eq_find_s2_region(points_sp, N);
    assert(r_idx == N, 'South pole should be in region N');
end

function test_boundaries
    N = 33;
    dim = 2;
    tol = eps * 2^5;
    s_regions = eq_regions(dim, N);
    [s_cap, ~] = eq_caps(dim, N);

    % 1. Test Cap Boundaries (Colatitudes)
    cap_points = zeros(2, length(s_cap));
    cap_points(2, :) = s_cap;

    % Points EXACTLY on s_cap[i] should be in the new cap (side="left" logic equivalent)
    % Expected regions for N=33: [1, 7, 16, 26, 32, 33]
    expected_cap_regions = [1, 7, 16, 26, 32, 33];
    r_idx_caps = eq_find_s2_region(cap_points, N);
    assert(isequal(r_idx_caps, expected_cap_regions), 'Cap boundary mismatch');

    % 2. Test Region Boundaries (Longitudes)
    % Using the start longitude of each region
    long_points = zeros(2, N);
    long_points(1, :) = squeeze(s_regions(1, 1, :));
    % Inside the colatitude band
    long_points(2, :) = squeeze(s_regions(2, 1, :)) + tol * 100;

    expected_long_regions = [1, 7, 2, 3, 4, 5, 6, 16, 8, 9, ...
                             10, 11, 12, 13, 14, 15, 26, 17, 18, 19, ...
                             20, 21, 22, 23, 24, 25, 32, 27, 28, 29, ...
                             30, 31, 33];
    r_idx_longs = eq_find_s2_region(long_points, N);
    assert(isequal(r_idx_longs, expected_long_regions), 'Longitude boundary mismatch');
end
