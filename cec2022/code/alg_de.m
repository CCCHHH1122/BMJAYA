function [best_f, best_x, curve] = alg_de(func_id, D, N, iters, lb, ub)
% ALG_DE  Differential Evolution: DE/rand/1/bin
% Signature matches your framework:
%   best_f = alg_de(func_id, D, N, iters, lb, ub)
%
% func_id : 1..12
% D       : dimension
% N       : population size
% iters   : MaxFEs
% lb, ub  : scalar bounds (e.g., -100, 100) or 1xD vectors

% ---- parameters (paper-friendly defaults) ----
F  = 0.5;
CR = 0.9;

Pop    = N;
MaxFEs = iters;

% ---- normalize bounds to 1xD row vectors ----
if isscalar(lb), lb = lb * ones(1, D); else, lb = reshape(lb, 1, []); end
if isscalar(ub), ub = ub * ones(1, D); else, ub = reshape(ub, 1, []); end

% ---- objective wrapper (CEC usually expects column vector) ----
fitfun = @(x) cec22_func(x', func_id);

% ---- init ----
pop = repmat(lb, Pop, 1) + rand(Pop, D) .* repmat(ub - lb, Pop, 1);
fit = zeros(Pop, 1);

fes = 0;
for i = 1:Pop
    fit(i) = fitfun(pop(i,:));
    fes = fes + 1;
end

[best_f, ib] = min(fit);
best_x = pop(ib,:);

% optional convergence curve (best-so-far)
curve = nan(MaxFEs, 1);
curve(1:fes) = best_f;

% ---- main loop ----
while fes < MaxFEs
    for i = 1:Pop
        if fes >= MaxFEs, break; end

        % pick r1,r2,r3 all distinct and != i
        idx = randperm(Pop, 3);
        while any(idx == i)
            idx = randperm(Pop, 3);
        end
        r1 = idx(1); r2 = idx(2); r3 = idx(3);

        % mutation rand/1
        v = pop(r1,:) + F * (pop(r2,:) - pop(r3,:));

        % boundary handling (clip)
        v = min(max(v, lb), ub);

        % binomial crossover
        u = pop(i,:);
        jrand = randi(D);
        for j = 1:D
            if rand < CR || j == jrand
                u(j) = v(j);
            end
        end

        fu = fitfun(u);
        fes = fes + 1;

        % greedy selection
        if fu <= fit(i)
            pop(i,:) = u;
            fit(i) = fu;

            if fu < best_f
                best_f = fu;
                best_x = u;
            end
        end

        curve(fes) = best_f;
    end
end

curve = curve(1:fes);
end
