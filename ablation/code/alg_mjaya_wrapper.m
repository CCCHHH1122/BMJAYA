function best_f = alg_mjaya_wrapper(func_id, D, N, iters, lb, ub)
K_req = 4;
MaxFEs = max(1, N + iters * N);
best_f = alg_mjaya_fixedN_budget(func_id, D, N, MaxFEs, lb, ub, K_req);
end
