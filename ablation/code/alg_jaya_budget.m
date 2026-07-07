function [best_f, curve] = alg_jaya_budget(func_id, D, N, MaxFEs, lb, ub)

lb = reshape(lb, 1, []);
ub = reshape(ub, 1, []);
if numel(lb) == 1, lb = repmat(lb, 1, D); end
if numel(ub) == 1, ub = repmat(ub, 1, D); end

curve = nan(MaxFEs, 1);

% ---- init ----
X = rand(N, D) .* (ub - lb) + lb;
fit = zeros(N, 1);
fes = 0;
best_so_far = inf;

for i = 1:N
    fit(i) = cec22_func(X(i, :), func_id);
    fes = fes + 1;
    best_so_far = min(best_so_far, fit(i));
    curve(fes) = best_so_far;

    if fes >= MaxFEs
        best_f = best_so_far;
        curve = fill_tail(curve);
        return;
    end
end

% ---- main loop ----
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

        Xj = X(i, :) + r1 .* (Xbest - X(i, :)) - r2 .* (Xworst - X(i, :));
        Xj = min(max(Xj, lb), ub);

        fj = cec22_func(Xj, func_id);
        fes = fes + 1;

        if fj < fit(i)
            X(i, :) = Xj;
            fit(i) = fj;
            best_so_far = min(best_so_far, fj);
        end

        curve(fes) = best_so_far;
    end
end

best_f = best_so_far;
curve = fill_tail(curve);
end

function curve = fill_tail(curve)
last = find(~isnan(curve), 1, 'last');
if isempty(last)
    return;
end
curve(last+1:end) = curve(last);
end
