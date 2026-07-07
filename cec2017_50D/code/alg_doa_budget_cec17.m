function [best_f, curve, fes] = alg_doa_budget_cec17(func_id, D, N, MaxFEs, lb, ub)
% Dhole Optimization Algorithm (DOA), paper-derived strict-budget adapter.
%
% Important:
% This implementation follows the equations and pseudocode in the DOA
% paper supplied by the user. Validate against the authors' source code
% before treating it as an exact official reproduction.

C1 = 1;
l = 25;
kappa = 0.5;
C3 = 3;

curve = nan(MaxFEs, 1);
fes = 0;

X = rand(N, D) .* (ub - lb) + lb;
fit = cec17_eval(X, func_id);

[curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fit, inf);

[~, bidx] = min(fit);
global_best_x = X(bidx, :);

T_equiv = max(1, ceil((MaxFEs - N) / N));
t = 0;

while fes < MaxFEs
    t = t + 1;

    [local_best_f, local_idx] = min(fit);
    local_best_x = X(local_idx, :);

    if local_best_f <= best_f
        global_best_x = local_best_x;
        best_f = local_best_f;
    end

    prey = (local_best_x + global_best_x) / 2;
    prey_f = max(cec17_eval(prey, func_id), eps);

    % Count the prey evaluation.
    [curve, fes, best_f] = cec17_append_batch_curve(curve, fes, prey_f, best_f);

    if fes >= MaxFEs
        break;
    end

    PMN = round(rand * 15 + 5);
    C2 = 1 - t / T_equiv;

    Xnew = X;

    for i = 1:N
        if fes >= MaxFEs
            break;
        end

        vocalization = rand;

        if vocalization < 0.5
            if PMN < 10
                % Eq. (6): search toward prey.
                Xcand = X(i, :) ...
                    + C2 .* rand(1, D) .* (prey - X(i, :));
            else
                % Eq. (8): encircling stage.
                z = randi(N);
                while N > 1 && z == i
                    z = randi(N);
                end

                Xcand = X(i, :) - X(z, :) + prey;
            end
        else
            % Eq. (4): hunting suitability.
            EF = rand;
            ps = C1 / (1 + exp(-kappa * (PMN - l))) * EF;

            % Eq. (10): prey size.
            S = C3 * rand * (fit(i) / prey_f);

            if S > (C3 + 1) / 2
                % Eq. (11) and Eq. (12): injury + repeated attack.
                Wprey = exp(-1 / max(S, eps)) .* local_best_x;

                Xcand = X(i, :) ...
                    + Wprey .* ps .* cos(2 * pi * rand(1, D)) ...
                    - ps .* sin(2 * pi * rand(1, D)) .* Wprey;
            else
                % Eq. (13): kill weak prey.
                Xcand = (X(i, :) - global_best_x) .* ps ...
                    + ps .* rand(1, D) .* X(i, :);
            end
        end

        Xcand = min(max(Xcand, lb), ub);

        fcand = cec17_eval(Xcand, func_id);
        [curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fcand, best_f);

        if fcand < fit(i)
            Xnew(i, :) = Xcand;
            fit(i) = fcand;

            if fcand <= best_f
                global_best_x = Xcand;
            end
        end
    end

    X = Xnew;
end

curve = cec17_fill_tail(curve);
end
