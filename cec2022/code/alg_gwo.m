function best_f = alg_gwo(func_id, D, N, iters, lb, ub)

X = rand(N,D).*(ub-lb) + lb;
fit = cec22_func(X, func_id);

T = max(iters, 1);

for t = 1:iters
    [~, idx] = sort(fit);
    alpha = X(idx(1),:);
    beta  = X(idx(2),:);
    delta = X(idx(3),:);

    a = 2 - 2*(t/T);

    r1 = rand(N,D); r2 = rand(N,D);
    A1 = 2*a.*r1 - a; C1 = 2.*r2;
    D1 = abs(C1.*alpha - X); X1 = alpha - A1.*D1;

    r1 = rand(N,D); r2 = rand(N,D);
    A2 = 2*a.*r1 - a; C2 = 2.*r2;
    D2 = abs(C2.*beta  - X); X2 = beta  - A2.*D2;

    r1 = rand(N,D); r2 = rand(N,D);
    A3 = 2*a.*r1 - a; C3 = 2.*r2;
    D3 = abs(C3.*delta - X); X3 = delta - A3.*D3;

    X = (X1 + X2 + X3) / 3;
    X = min(max(X, lb), ub);

    fit = cec22_func(X, func_id);
end

best_f = min(fit);
end
