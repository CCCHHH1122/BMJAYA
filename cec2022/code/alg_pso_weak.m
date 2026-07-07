function best_f = alg_pso_weak(func_id, D, N, iters, lb, ub)

c1 = 2.0; c2 = 2.0;

X = rand(N,D).*(ub-lb) + lb;
V = zeros(N,D);

fit = cec22_func(X, func_id);
P   = X;
Pfit= fit;

[best_f, gidx] = min(fit);
G = X(gidx,:);

T = max(iters, 1);

for t = 1:iters
    w = 0.9 - 0.6*(t/T);          % linearly decreasing inertia
    r1 = rand(N,D);
    r2 = rand(N,D);

    V = w.*V + c1.*r1.*(P - X) + c2.*r2.*(G - X);
    X = X + V;

    % simple clip boundary (weaker, common baseline)
    X = min(max(X, lb), ub);

    fit = cec22_func(X, func_id);

    improved = fit < Pfit;
    P(improved,:) = X(improved,:);
    Pfit(improved) = fit(improved);

    [cur_best, gidx] = min(fit);
    if cur_best < best_f
        best_f = cur_best;
        G = X(gidx,:);
    end
end

end
