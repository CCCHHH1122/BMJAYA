function [best_f, curve] = alg_bjaya_budget(func_id, D, N, MaxFEs, lb, ub)
% BMRJAYA under strict MaxFEs budget (single population)
% - Stage I: JAYA
% - Stage II: BMR repair triggered ONLY if JAYA fails
% - strict MaxFEs: never exceed evaluation budget
% - outputs convergence curve by FEs: curve(t) = best-so-far after t-th evaluation
%
% Signature:
%   [best_f, curve] = alg_bmrjaya_budget_curve(func_id, D, N, MaxFEs, lb, ub)

% ---- bounds normalize ----
lb = reshape(lb,1,[]);
ub = reshape(ub,1,[]);
if numel(lb)==1, lb = repmat(lb,1,D); end
if numel(ub)==1, ub = repmat(ub,1,D); end

curve = nan(MaxFEs,1);

% ----- init -----
X = rand(N,D).*(ub-lb) + lb;

% 如果你的 cec22_func 支持矩阵输入，这句可用：fit = cec22_func(X, func_id);
% 为了兼容向量版，这里逐行评估：
fit = zeros(N,1);
fes = 0;
best_so_far = inf;

for i = 1:N
    fit(i) = cec22_func(X(i,:), func_id);
    fes = fes + 1;
    best_so_far = min(best_so_far, fit(i));
    curve(fes) = best_so_far;
    if fes >= MaxFEs
        best_f = best_so_far;
        curve = fill_tail(curve);
        return;
    end
end

% ----- main loop -----
while fes < MaxFEs

    % population stats
    [~, bidx] = min(fit);
    [~, widx] = max(fit);
    Xbest  = X(bidx,:);
    Xworst = X(widx,:);
    Xmean  = mean(X,1);

    for i = 1:N

        % ---- Stage I: JAYA ----
        if fes >= MaxFEs
            best_f = best_so_far;
            curve = fill_tail(curve);
            return;
        end

        r1 = rand(1,D);
        r2 = rand(1,D);

        Xj = X(i,:) + r1.*(Xbest - X(i,:)) - r2.*(Xworst - X(i,:));
        Xj = min(max(Xj, lb), ub);

        fj = cec22_func(Xj, func_id);
        fes = fes + 1;

        if fj < fit(i)
            X(i,:) = Xj;
            fit(i) = fj;
            best_so_far = min(best_so_far, fj);
            curve(fes) = best_so_far;
            continue; % skip BMR
        else
            % JAYA没改进，也要记录best_so_far（保持按FEs长度对齐）
            curve(fes) = best_so_far;
        end

        % ---- Stage II: BMR (triggered only if JAYA fails) ----
        if fes >= MaxFEs
            best_f = best_so_far;
            curve = fill_tail(curve);
            return;
        end

        ridx  = randi(N);  % allow self-selection (consistent with your BMJAYA)
        Xrand = X(ridx,:);

        Xb = Xbest + (Xmean - X(i,:)) + (Xrand - X(i,:));
        Xb = min(max(Xb, lb), ub);

        fb = cec22_func(Xb, func_id);
        fes = fes + 1;

        if fb < fit(i)
            X(i,:) = Xb;
            fit(i) = fb;
            best_so_far = min(best_so_far, fb);
        end
        curve(fes) = best_so_far;
    end
end

best_f = best_so_far;
curve = fill_tail(curve);
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
