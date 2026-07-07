function best_f = alg_bmjaya_fixedN_budget(func_id, D, N, MaxFEs, lb, ub, K_req)
% BMJAYA (fixed total population N) under strict MaxFEs budget
% - multi-subpopulation sizes sum to N
% - per-individual: JAYA first; if improved, skip BMR; else BMR
% - elite migration EVERY generation
% - strict MaxFEs budget
%
% Signature:
%   best_f = alg_bmjaya_fixedN_budget(func_id, D, N, MaxFEs, lb, ub, K_req)
%
% NOTE:
%   Requires: cec22_func(x, func_id) in your path (either .m or mex).

lb = reshape(lb,1,[]);
ub = reshape(ub,1,[]);

% If scalar bounds, expand to 1xD (recommended for safety)
if numel(lb) == 1, lb = repmat(lb, 1, D); end
if numel(ub) == 1, ub = repmat(ub, 1, D); end

K = min(max(1, K_req), N);
sizes = split_sizes(N, K);

% -------- init subpops --------
subX = cell(K,1);
subF = cell(K,1);
fes  = 0;

for k = 1:K
    nk = sizes(k);
    X = rand(nk,D).*(ub-lb) + lb;
    f = zeros(nk,1);
    for i = 1:nk
        f(i) = cec22_func(X(i,:), func_id);
    end
    fes = fes + nk;

    subX{k} = X;
    subF{k} = f;
end

% -------- init global best --------
gbest_f = inf;
gbest_x = [];

for k = 1:K
    [bk, bi] = min(subF{k});
    if bk < gbest_f
        gbest_f = bk;
        gbest_x = subX{k}(bi,:);
    end
end

% ================= main loop =================
while fes < MaxFEs

    % ---- local updates (one generation) ----
    for k = 1:K
        [subX{k}, subF{k}, fes, gbest_x, gbest_f] = ...
            bmr_jaya_step_budget_sync_gbest( ...
                subX{k}, subF{k}, func_id, ...
                lb, ub, fes, MaxFEs, gbest_x, gbest_f);

        if fes >= MaxFEs
            best_f = gbest_f;
            return;
        end
    end

    % ---- elite migration (same as Python) ----
    for k = 1:K
        [wk, wi] = max(subF{k});
        if gbest_f < wk
            subX{k}(wi,:) = gbest_x;
            subF{k}(wi)   = gbest_f;
        end
    end
end

best_f = gbest_f;
end


% =========================================================
% One-generation update for ONE subpopulation
% - JAYA first; if improved => accept + skip BMR
% - else try BMR
% - strict MaxFEs check
% - IMMEDIATE global-best sync (critical for budget stopping)
% =========================================================
function [X, fit, fes, gbest_x, gbest_f] = ...
    bmr_jaya_step_budget_sync_gbest(X, fit, func_id, lb, ub, fes, MaxFEs, gbest_x, gbest_f)

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

        % immediate gbest sync
        if fj < gbest_f
            gbest_f = fj;
            gbest_x = Xj;
        end
        continue; % skip BMR
    end

    % ---- Stage 2: BMR (only if JAYA fails) ----
    if fes >= MaxFEs; return; end

    ridx  = randi(nk);     % (Python allows self-selection; keep consistent)
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
