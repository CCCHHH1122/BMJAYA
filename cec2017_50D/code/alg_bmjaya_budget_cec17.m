function [best_f, curve, fes] = alg_bmjaya_budget_cec17(func_id, D, N, MaxFEs, lb, ub)
% BMJAYA with strict MaxFEs.
%
% Preserves the user's original BMJAYA rule:
%   - fixed total population N
%   - K = 4 subpopulations
%   - JAYA first
%   - conditional BMR repair only when JAYA fails
%   - BMR: Xbest + (Xmean-Xi) + (Xrand-Xi)
%   - self-selection is allowed for Xrand
%   - elite sharing every generation

K_req = 4;

curve = nan(MaxFEs, 1);
fes = 0;

K = min(max(1, K_req), N);
sizes = split_sizes(N, K);

subX = cell(K, 1);
subF = cell(K, 1);

gbest_f = inf;
gbest_x = zeros(1, D);

for k = 1:K
    nk = sizes(k);

    X = rand(nk, D) .* (ub - lb) + lb;
    f = cec17_eval(X, func_id);

    [curve, fes, gbest_f] = cec17_append_batch_curve(curve, fes, f, gbest_f);

    [local_best, local_idx] = min(f);
    if local_best <= gbest_f
        gbest_x = X(local_idx, :);
        gbest_f = local_best;
    end

    subX{k} = X;
    subF{k} = f;
end

while fes < MaxFEs
    for k = 1:K
        [subX{k}, subF{k}, fes, gbest_x, gbest_f, curve] = ...
            local_step(subX{k}, subF{k}, func_id, lb, ub, ...
                fes, MaxFEs, gbest_x, gbest_f, curve);

        if fes >= MaxFEs
            best_f = gbest_f;
            curve = cec17_fill_tail(curve);
            return;
        end
    end

    % Elite sharing every generation.
    for k = 1:K
        [worst_f, worst_idx] = max(subF{k});

        if gbest_f < worst_f
            subX{k}(worst_idx, :) = gbest_x;
            subF{k}(worst_idx) = gbest_f;
        end
    end
end

best_f = gbest_f;
curve = cec17_fill_tail(curve);
end

function [X, fit, fes, gbest_x, gbest_f, curve] = ...
    local_step(X, fit, func_id, lb, ub, fes, MaxFEs, gbest_x, gbest_f, curve)

[nk, D] = size(X);

Xmean = mean(X, 1);

[~, bidx] = min(fit);
[~, widx] = max(fit);

Xbest = X(bidx, :);
Xworst = X(widx, :);

for i = 1:nk
    if fes >= MaxFEs
        return;
    end

    r1 = rand(1, D);
    r2 = rand(1, D);

    Xj = X(i, :) ...
       + r1 .* (Xbest - X(i, :)) ...
       - r2 .* (Xworst - X(i, :));

    Xj = min(max(Xj, lb), ub);

    fj = cec17_eval(Xj, func_id);
    [curve, fes, gbest_f] = cec17_append_batch_curve(curve, fes, fj, gbest_f);

    if fj < fit(i)
        X(i, :) = Xj;
        fit(i) = fj;

        if fj <= gbest_f
            gbest_f = fj;
            gbest_x = Xj;
        end

        continue;
    end

    if fes >= MaxFEs
        return;
    end

    % Original conditional BMR repair.
    ridx = randi(nk);  % self-selection intentionally allowed
    Xrand = X(ridx, :);

    Xb = Xbest ...
       + (Xmean - X(i, :)) ...
       + (Xrand - X(i, :));

    Xb = min(max(Xb, lb), ub);

    fb = cec17_eval(Xb, func_id);
    [curve, fes, gbest_f] = cec17_append_batch_curve(curve, fes, fb, gbest_f);

    if fb < fit(i)
        X(i, :) = Xb;
        fit(i) = fb;

        if fb <= gbest_f
            gbest_f = fb;
            gbest_x = Xb;
        end
    end
end
end

function sizes = split_sizes(N, K)
base = floor(N / K);
remainder = N - base * K;

sizes = base * ones(K, 1);
sizes(1:remainder) = sizes(1:remainder) + 1;
end
