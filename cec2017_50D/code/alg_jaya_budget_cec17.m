function [best_f, curve, fes] = alg_jaya_budget_cec17(func_id, D, N, MaxFEs, lb, ub)
% Standard JAYA with strict MaxFEs.
curve = nan(MaxFEs, 1);
fes = 0;

X = rand(N, D) .* (ub - lb) + lb;
fit = cec17_eval(X, func_id);
[curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fit, inf);

while fes < MaxFEs
    [~, bidx] = min(fit);
    [~, widx] = max(fit);

    Xbest = X(bidx, :);
    Xworst = X(widx, :);

    for i = 1:N
        if fes >= MaxFEs
            break;
        end

        r1 = rand(1, D);
        r2 = rand(1, D);

        Xnew = X(i, :) ...
            + r1 .* (Xbest - X(i, :)) ...
            - r2 .* (Xworst - X(i, :));

        Xnew = min(max(Xnew, lb), ub);

        fnew = cec17_eval(Xnew, func_id);
        [curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fnew, best_f);

        if fnew < fit(i)
            X(i, :) = Xnew;
            fit(i) = fnew;
        end
    end
end

curve = cec17_fill_tail(curve);
end
