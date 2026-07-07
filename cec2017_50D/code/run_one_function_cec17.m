function [results, sampled_curves, fes_used] = run_one_function_cec17( ...
    algo_fn, func_id, D, N, MaxFEs, Runs, lb, ub, fes_show, baseSeed, alg_id)
% Execute one optimizer on one CEC2017 function over independent runs.

results = nan(Runs, 1);
sampled_curves = nan(Runs, numel(fes_show), 'single');
fes_used = zeros(Runs, 1);

for run_id = 1:Runs
    seed = baseSeed + D * 1000000 + func_id * 10000 + alg_id * 100 + run_id;
    rng(seed, 'twister');

    [best_f, curve, fes] = algo_fn(func_id, D, N, MaxFEs, lb, ub);

    assert(fes == MaxFEs, ...
        'F%d: optimizer used %d FEs instead of %d.', ...
        func_id, fes, MaxFEs);

    assert(all(isfinite(curve)), ...
        'F%d: optimizer returned NaN or Inf in convergence curve.', func_id);

    assert(all(diff(curve) <= 1e-10), ...
        'F%d: convergence curve is not monotonically non-increasing.', func_id);

    results(run_id) = best_f;
    sampled_curves(run_id, :) = single(curve(fes_show));
    fes_used(run_id) = fes;
end
end
