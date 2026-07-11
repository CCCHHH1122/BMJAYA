function [best_f, best_x, curve] = BMJAYA(objfun, D, N, MaxFEs, lb, ub, K_req)
% BMJAYA  Parameter-free multi-strategy enhanced JAYA optimizer.
%
%   [best_f, best_x, curve] = BMJAYA(objfun, D, N, MaxFEs, lb, ub)
%   minimizes the objective function handle objfun over a D-dimensional
%   bounded search space.
%
%   Inputs:
%     objfun  - function handle, for example @(x) sum(x.^2)
%     D       - problem dimension
%     N       - total population size
%     MaxFEs  - maximum number of function evaluations
%     lb      - lower bound, scalar or 1-by-D vector
%     ub      - upper bound, scalar or 1-by-D vector
%     K_req   - optional number of subpopulations; default is 4
%
%   Outputs:
%     best_f  - best objective value found
%     best_x  - best solution vector found
%     curve   - best-so-far value after each function evaluation
%
%   Core BMJAYA mechanisms:
%     1. fixed total population split into multiple subpopulations;
%     2. JAYA update attempted first for each individual;
%     3. conditional Best-Mean-Random (BMR) repair only if JAYA fails;
%     4. global elite sharing after each generation.

if nargin < 7 || isempty(K_req)
    K_req = 4;
end

lb = reshape(lb, 1, []);
ub = reshape(ub, 1, []);

if numel(lb) == 1
    lb = repmat(lb, 1, D);
end
if numel(ub) == 1
    ub = repmat(ub, 1, D);
end

if numel(lb) ~= D || numel(ub) ~= D
    error('Bounds must be scalars or vectors with D elements.');
end
if N < 1 || MaxFEs < 1
    error('N and MaxFEs must be positive integers.');
end

K = min(max(1, round(K_req)), N);
sizes = split_sizes(N, K);

subX = cell(K, 1);
subF = cell(K, 1);
curve = nan(MaxFEs, 1);

fes = 0;
best_f = inf;
best_x = zeros(1, D);

for k = 1:K
    nk = sizes(k);
    X = rand(nk, D) .* (ub - lb) + lb;
    f = inf(nk, 1);

    for i = 1:nk
        if fes >= MaxFEs
            curve = fill_curve_tail(curve);
            return;
        end

        f(i) = objfun(X(i, :));
        fes = fes + 1;

        if f(i) < best_f
            best_f = f(i);
            best_x = X(i, :);
        end
        curve(fes) = best_f;
    end

    subX{k} = X;
    subF{k} = f;
end

while fes < MaxFEs
    for k = 1:K
        [subX{k}, subF{k}, fes, best_x, best_f, curve] = ...
            local_step(subX{k}, subF{k}, objfun, lb, ub, ...
            fes, MaxFEs, best_x, best_f, curve);

        if fes >= MaxFEs
            curve = fill_curve_tail(curve);
            return;
        end
    end

    for k = 1:K
        [worst_f, worst_idx] = max(subF{k});
        if best_f < worst_f
            subX{k}(worst_idx, :) = best_x;
            subF{k}(worst_idx) = best_f;
        end
    end
end

curve = fill_curve_tail(curve);
end

function [X, fit, fes, best_x, best_f, curve] = ...
    local_step(X, fit, objfun, lb, ub, fes, MaxFEs, best_x, best_f, curve)

[nk, D] = size(X);
Xmean = mean(X, 1);

[~, best_idx] = min(fit);
[~, worst_idx] = max(fit);
Xbest = X(best_idx, :);
Xworst = X(worst_idx, :);

for i = 1:nk
    if fes >= MaxFEs
        return;
    end

    r1 = rand(1, D);
    r2 = rand(1, D);

    Xj = X(i, :) ...
        + r1 .* (Xbest - X(i, :)) ...
        - r2 .* (Xworst - X(i, :));
    Xj = bound_clip(Xj, lb, ub);

    fj = objfun(Xj);
    fes = fes + 1;

    if fj < best_f
        best_f = fj;
        best_x = Xj;
    end
    curve(fes) = best_f;

    if fj < fit(i)
        X(i, :) = Xj;
        fit(i) = fj;
        continue;
    end

    if fes >= MaxFEs
        return;
    end

    random_idx = randi(nk);
    Xrand = X(random_idx, :);

    Xb = Xbest + (Xmean - X(i, :)) + (Xrand - X(i, :));
    Xb = bound_clip(Xb, lb, ub);

    fb = objfun(Xb);
    fes = fes + 1;

    if fb < best_f
        best_f = fb;
        best_x = Xb;
    end
    curve(fes) = best_f;

    if fb < fit(i)
        X(i, :) = Xb;
        fit(i) = fb;
    end
end
end

function x = bound_clip(x, lb, ub)
x = min(max(x, lb), ub);
end

function sizes = split_sizes(N, K)
base = floor(N / K);
remainder = N - base * K;
sizes = base * ones(K, 1);
sizes(1:remainder) = sizes(1:remainder) + 1;
end

function curve = fill_curve_tail(curve)
last = find(~isnan(curve), 1, 'last');
if ~isempty(last) && last < numel(curve)
    curve(last+1:end) = curve(last);
end
end
