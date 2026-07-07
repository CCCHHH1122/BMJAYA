function best_f = alg_bmjaya_wrapper(func_id, D, N, iters, lb, ub)
% Wrapper for BMJAYA to keep iter-based interface:
%   best_f = algo_fn(func_id, D, N, iters, lb, ub)
%
% Compatible with your iters_from_maxfes(MaxFEs, N):
%   iters = floor((MaxFEs - N) / N)
%
% So inverse mapping (consistent budget) is:
%   MaxFEs = N + iters * N

K_req = 4; % fixed internal subpopulation setting used in the manuscript

% Inverse of iters_from_maxfes (budget aligned with other algorithms)
MaxFEs = max(1, N + iters * N);

best_f = alg_bmjaya_fixedN_budget(func_id, D, N, MaxFEs, lb, ub, K_req);
end
