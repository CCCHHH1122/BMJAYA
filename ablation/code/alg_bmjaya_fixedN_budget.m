function [best_f, curve] = alg_bmjaya_fixedN_budget(func_id, D, N, MaxFEs, lb, ub, K_req)
% BMJAYA (fixed total population N) under strict MaxFEs budget
% Mechanisms:
%   - Multipopulation Cooperative Search (K subpops, sizes sum to N)
%   - Elite Sharing (global best migration every generation)
%   - BMR mechanism (triggered if JAYA fails for an individual)
% Budget:
%   - strict MaxFEs: never exceed evaluation budget
% Output:
%   - curve(t) = global best-so-far after t-th function evaluation (FEs)
%
% Signature:
%   [best_f, curve] = alg_bmjaya_fixedN_budget_curve(func_id, D, N, MaxFEs, lb, ub, K_req)

lb = reshape(lb,1,[]);
ub = reshape(ub,1,[]);
if numel(lb) == 1, lb = repmat(lb, 1, D); end
if numel(ub) == 1, ub = repmat(ub, 1, D); end

curve = nan(MaxFEs,1);

K = min(max(1, K_req), N);
sizes = split_sizes(N, K);

% -------- init subpops --------
subX = cell(K,1);
subF = cell(K,1);
fes  = 0;

gbest_f = inf;
gbest_x = zeros(1,D);

for k = 1:K
    nk = sizes(k);
    X = rand(nk,D).*(ub-lb) + lb;
    f = zeros(nk,1);

    for i = 1:nk
        f(i) = cec22_func(X(i,:), func_id);
        fes = fes + 1;

        if f(i) < gbest_f
            gbest_f = f(i);
            gbest_x = X(i,:);
        end

        if fes <= MaxFEs
            curve(fes) = gbest_f;
        end

        if fes >= MaxFEs
            subX{k} = X;
            subF{k} = f;
            best_f = gbest_f;
            curve = fill_tail(curve);
            return;
        end
    end

    subX{k} = X;
    subF{k} = f;
end

% ================= main loop =================
while fes < MaxFEs

    % ---- local updates (one generation per subpopulation) ----
    for k = 1:K
        [subX{k}, subF{k}, fes, gbest_x, gbest_f, curve] = ...
            bmr_jaya_step_budget_sync_gbest_curve( ...
                subX{k}, subF{k}, func_id, ...
                lb, ub, fes, MaxFEs, gbest_x, gbest_f, curve);

        if fes >= MaxFEs
            best_f = gbest_f;
            curve = fill_tail(curve);
            return;
        end
    end

    % ---- elite migration (global elite sharing) ----
    for k = 1:K
        [wk, wi] = max(subF{k});
        if gbest_f < wk
            subX{k}(wi,:) = gbest_x;
            subF{k}(wi)   = gbest_f;
        end
    end
end

best_f = gbest_f;
curve = fill_tail(curve);
end


% =========================================================
% One-generation update for ONE subpopulation
% - JAYA first; if improved => accept + skip BMR
% - else try BMR
% - strict MaxFEs check
% - IMMEDIATE global-best sync + curve logging
% =========================================================
function [X, fit, fes, gbest_x, gbest_f, curve] = ...
    bmr_jaya_step_budget_sync_gbest_curve(X, fit, func_id, lb, ub, fes, MaxFEs, gbest_x, gbest_f, curve)

[nk, D] = size(X);

Xmean = mean(X, 1);

[~, bidx] = min(fit);
[~, widx] = max(fit);
Xbest  = X(bidx,:);
Xworst = X(widx,:);

for i = 1:nk

    % ---- Stage 1: JAYA ----
    if fes >= MaxFEs; return; end

    r1 = rand(1, D);
    r2 = rand(1, D);

    Xj = X(i,:) + r1.*(Xbest - X(i,:)) - r2.*(Xworst - X(i,:));
    Xj = min(max(Xj, lb), ub);

    fj = cec22_func(Xj, func_id);
    fes = fes + 1;

    if fj < fit(i)
        X(i,:) = Xj;
        fit(i) = fj;

        if fj < gbest_f
            gbest_f = fj;
            gbest_x = Xj;
        end
        curve(fes) = gbest_f;
        continue; % skip BMR
    else
        curve(fes) = gbest_f;
    end

    % ---- Stage 2: BMR (only if JAYA fails) ----
    if fes >= MaxFEs; return; end

    ridx  = randi(nk);     % allow self-selection (consistent with your BMJAYA)
    Xrand = X(ridx,:);

    Xb = Xbest + (Xmean - X(i,:)) + (Xrand - X(i,:));
    Xb = min(max(Xb, lb), ub);

    fb = cec22_func(Xb, func_id);
    fes = fes + 1;

    if fb < fit(i)
        X(i,:) = Xb;
        fit(i) = fb;

        if fb < gbest_f
            gbest_f = fb;
            gbest_x = Xb;
        end
    end
    curve(fes) = gbest_f;
end
end


% =========================================================
% Split N into K near-equal sizes summing to N
% =========================================================
function sizes = split_sizes(N, K)
base = floor(N / K);
remd = N - base*K;
sizes = base * ones(K, 1);
for i = 1:remd
    sizes(i) = sizes(i) + 1;
end
end

% ---- helper: fill tail NaNs if any ----
function curve = fill_tail(curve)
last = find(~isnan(curve), 1, 'last');
if isempty(last)
    curve(:) = inf;
else
    curve(last+1:end) = curve(last);
end
end
