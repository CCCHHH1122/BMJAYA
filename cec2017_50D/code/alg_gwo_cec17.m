function [best_f, curve, fes] = alg_gwo_cec17(func_id, D, N, MaxFEs, lb, ub)
% Standard GWO baseline with strict MaxFEs.
curve = nan(MaxFEs, 1);
fes = 0;

X = rand(N, D) .* (ub - lb) + lb;
fit = cec17_eval(X, func_id);
[curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fit, inf);

max_steps = max(1, ceil((MaxFEs - fes) / N));
step = 0;

while fes < MaxFEs
    step = step + 1;

    [~, idx] = sort(fit);
    alpha = X(idx(1), :);
    beta  = X(idx(2), :);
    delta = X(idx(3), :);

    a = 2 - 2 * step / max_steps;

    r1 = rand(N, D); r2 = rand(N, D);
    A1 = 2 * a .* r1 - a; C1 = 2 .* r2;
    X1 = alpha - A1 .* abs(C1 .* alpha - X);

    r1 = rand(N, D); r2 = rand(N, D);
    A2 = 2 * a .* r1 - a; C2 = 2 .* r2;
    X2 = beta - A2 .* abs(C2 .* beta - X);

    r1 = rand(N, D); r2 = rand(N, D);
    A3 = 2 * a .* r1 - a; C3 = 2 .* r2;
    X3 = delta - A3 .* abs(C3 .* delta - X);

    Xcand = min(max((X1 + X2 + X3) / 3, lb), ub);

    nEval = min(N, MaxFEs - fes);
    fnew = cec17_eval(Xcand(1:nEval, :), func_id);
    [curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fnew, best_f);

    X(1:nEval, :) = Xcand(1:nEval, :);
    fit(1:nEval) = fnew;
end

curve = cec17_fill_tail(curve);
end
