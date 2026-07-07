function results = run_one_function(algo_fn, func_id, D, N, MaxFEs, Runs, lb, ub)

iters = iters_from_maxfes(MaxFEs, N);
call_arg = budget_argument_for_algorithm(algo_fn, MaxFEs, iters);
results = zeros(Runs, 1);

for r = 1:Runs
    rng(1 + r, 'twister'); % reproducible
    best_f = algo_fn(func_id, D, N, call_arg, lb, ub);
    results(r) = best_f;
end

end

function call_arg = budget_argument_for_algorithm(algo_fn, MaxFEs, iters)
% Wrapper-based variants receive iteration counts and restore MaxFEs inside.
% Budget-based implementations already use the fourth argument as MaxFEs.
maxfes_algorithms = { ...
    'alg_jaya_budget', ...
    'alg_bjaya_budget' ...
    };

if ismember(func2str(algo_fn), maxfes_algorithms)
    call_arg = MaxFEs;
else
    call_arg = iters;
end
end
