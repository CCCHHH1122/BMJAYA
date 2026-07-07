function [best_f, curve] = alg_mjaya_fixedN_budget(func_id, D, N, MaxFEs, lb, ub, K_req)

lb = reshape(lb,1,[]);
ub = reshape(ub,1,[]);
if numel(lb) == 1, lb = repmat(lb, 1, D); end
if numel(ub) == 1, ub = repmat(ub, 1, D); end

curve = nan(MaxFEs,1);

K = min(max(1, K_req), N);
sizes = split_sizes(N, K);

subX = cell(K,1);
subF = cell(K,1);
fes  = 0;

for k = 1:K
    nk = sizes(k);
    X = rand(nk,D).*(ub-lb) + lb;

    f = zeros(nk,1);
    for i = 1:nk
        f(i) = cec22_func(X(i,:), func_id);
        fes = fes + 1;
        % 暂时先不写curve，等全局best初始化完再补
    end
    subX{k} = X;
    subF{k} = f;
end

% init global best
gbest_f = inf;
gbest_x = zeros(1,D);
for k = 1:K
    [bk, bi] = min(subF{k});
    if bk < gbest_f
        gbest_f = bk;
        gbest_x = subX{k}(bi,:);
    end
end
curve(1:fes) = gbest_f;

while fes < MaxFEs

    for k = 1:K
        [subX{k}, subF{k}, fes, gbest_x, gbest_f, curve] = ...
            jaya_step_budget_sync_gbest_curve( ...
                subX{k}, subF{k}, func_id, lb, ub, fes, MaxFEs, gbest_x, gbest_f, curve);

        if fes >= MaxFEs
            best_f = gbest_f;
            curve = fill_tail(curve);
            return;
        end
    end

    % elite sharing
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


function [X, fit, fes, gbest_x, gbest_f, curve] = ...
    jaya_step_budget_sync_gbest_curve(X, fit, func_id, lb, ub, fes, MaxFEs, gbest_x, gbest_f, curve)

[nk, D] = size(X);
[~, bidx] = min(fit);
[~, widx] = max(fit);
Xbest  = X(bidx,:);
Xworst = X(widx,:);

for i = 1:nk
    if fes >= MaxFEs, return; end

    r1 = rand(1,D);
    r2 = rand(1,D);

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
    end

    curve(fes) = gbest_f;
end
end

function sizes = split_sizes(N, K)
base = floor(N / K);
remd = N - base*K;
sizes = base * ones(K, 1);
for i = 1:remd
    sizes(i) = sizes(i) + 1;
end
end

function curve = fill_tail(curve)
last = find(~isnan(curve),1,'last');
curve(last+1:end) = curve(last);
end
