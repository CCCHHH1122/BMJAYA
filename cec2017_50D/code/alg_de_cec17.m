function [best_f, curve, fes] = alg_de_cec17(func_id, D, N, MaxFEs, lb, ub)
% DE/rand/1/bin with strict MaxFEs.
F = 0.5;
CR = 0.9;

curve = nan(MaxFEs, 1);
fes = 0;

pop = rand(N, D) .* (ub - lb) + lb;
fit = cec17_eval(pop, func_id);
[curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fit, inf);

while fes < MaxFEs
    for i = 1:N
        if fes >= MaxFEs
            break;
        end

        idx = randperm(N, 3);
        while any(idx == i)
            idx = randperm(N, 3);
        end

        v = pop(idx(1), :) + F * (pop(idx(2), :) - pop(idx(3), :));
        v = min(max(v, lb), ub);

        u = pop(i, :);
        jrand = randi(D);

        for j = 1:D
            if rand < CR || j == jrand
                u(j) = v(j);
            end
        end

        fu = cec17_eval(u, func_id);
        [curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fu, best_f);

        if fu <= fit(i)
            pop(i, :) = u;
            fit(i) = fu;
        end
    end
end

curve = cec17_fill_tail(curve);
end
