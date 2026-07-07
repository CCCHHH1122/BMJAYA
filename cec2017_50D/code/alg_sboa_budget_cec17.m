function [best_f, curve, fes] = alg_sboa_budget_cec17(func_id, D, N, MaxFEs, lb, ub)
% Secretary Bird Optimization Algorithm (SBOA) with strict MaxFEs.
curve = nan(MaxFEs, 1);
fes = 0;

X = rand(N, D) .* (ub - lb) + lb;
fit = cec17_eval(X, func_id);

[curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fit, inf);
[~, bidx] = min(fit);
best_x = X(bidx, :);

T_equiv = max(1, ceil((MaxFEs - N) / (2 * N)));
t = 0;

while fes < MaxFEs
    t = t + 1;
    ratio = min(t / T_equiv, 1);

    % Phase P1: hunting behavior
    for i = 1:N
        if fes >= MaxFEs
            break;
        end

        if t < T_equiv / 3
            idx = randperm(N, 2);
            Xnew = X(i, :) ...
                + (X(idx(1), :) - X(idx(2), :)) .* rand(1, D);

        elseif t < 2 * T_equiv / 3
            RB = randn(1, D);

            Xnew = best_x ...
                + exp(ratio^4) .* (RB - 0.5) .* (best_x - X(i, :));

        else
            RL = 0.5 .* levy_flight(D);

            Xnew = best_x ...
                + ((1 - ratio)^(2 * ratio)) .* X(i, :) .* RL;
        end

        Xnew = min(max(Xnew, lb), ub);

        fnew = cec17_eval(Xnew, func_id);
        [curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fnew, best_f);

        if fnew < fit(i)
            X(i, :) = Xnew;
            fit(i) = fnew;

            if fnew <= best_f
                best_x = Xnew;
            end
        end
    end

    % Phase P2: escaping from predators
    for i = 1:N
        if fes >= MaxFEs
            break;
        end

        if rand < 0.5
            RB = randn(1, D);

            Xnew = best_x ...
                + (2 .* RB - 1) .* ((1 - ratio)^2) .* X(i, :);

        else
            K = round(1 + rand);
            R2 = randn(1, D);

            random_idx = randi(N);
            Xrandom = X(random_idx, :);

            Xnew = X(i, :) + R2 .* (Xrandom - K .* X(i, :));
        end

        Xnew = min(max(Xnew, lb), ub);

        fnew = cec17_eval(Xnew, func_id);
        [curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fnew, best_f);

        if fnew < fit(i)
            X(i, :) = Xnew;
            fit(i) = fnew;

            if fnew <= best_f
                best_x = Xnew;
            end
        end
    end
end

curve = cec17_fill_tail(curve);
end

function step = levy_flight(D)
beta = 1.5;
scale = 0.01;

sigma = (gamma(1 + beta) * sin(pi * beta / 2) / ...
    (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);

u = randn(1, D) .* sigma;
v = randn(1, D);

step = scale .* u ./ (abs(v).^(1 / beta) + eps);
end
