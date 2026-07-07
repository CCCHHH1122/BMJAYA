function [best_f, curve, fes] = alg_pso_weak_cec17(func_id, D, N, MaxFEs, lb, ub)
% Standard PSO baseline with strict MaxFEs.
c1 = 2.0;
c2 = 2.0;

curve = nan(MaxFEs, 1);
fes = 0;

X = rand(N, D) .* (ub - lb) + lb;
V = zeros(N, D);

fit = cec17_eval(X, func_id);
[curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fit, inf);

P = X;
Pfit = fit;
[~, gidx] = min(Pfit);
G = P(gidx, :);

max_steps = max(1, ceil((MaxFEs - fes) / N));
step = 0;

while fes < MaxFEs
    step = step + 1;
    w = 0.9 - 0.6 * step / max_steps;

    r1 = rand(N, D);
    r2 = rand(N, D);

    Vcand = w .* V + c1 .* r1 .* (P - X) + c2 .* r2 .* (G - X);
    Xcand = min(max(X + Vcand, lb), ub);

    nEval = min(N, MaxFEs - fes);
    fnew = cec17_eval(Xcand(1:nEval, :), func_id);
    [curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fnew, best_f);

    X(1:nEval, :) = Xcand(1:nEval, :);
    V(1:nEval, :) = Vcand(1:nEval, :);

    improved = fnew < Pfit(1:nEval);
    idx = find(improved);

    if ~isempty(idx)
        P(idx, :) = X(idx, :);
        Pfit(idx) = fnew(improved);
    end

    [pbest, gidx] = min(Pfit);
    G = P(gidx, :);
    best_f = min(best_f, pbest);
end

curve = cec17_fill_tail(curve);
end
