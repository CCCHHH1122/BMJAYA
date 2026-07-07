function best_f = alg_jaya_budget(func_id, D, N, MaxFEs, lb, ub)
% JAYA under strict MaxFEs budget (single population)
% - strict MaxFEs: never exceed evaluation budget
%
% Signature:
%   best_f = alg_jaya_budget(func_id, D, N, MaxFEs, lb, ub)

% ---- bounds normalize ----
lb = reshape(lb,1,[]);
ub = reshape(ub,1,[]);
if numel(lb)==1, lb = repmat(lb,1,D); end
if numel(ub)==1, ub = repmat(ub,1,D); end

% ---- init ----
X = rand(N,D).*(ub-lb) + lb;
fit = cec22_func(X, func_id);   % N evaluations
fes = N;

% ---- main loop ----
while fes < MaxFEs

    % population statistics (current generation)
    [~, bidx] = min(fit);
    [~, widx] = max(fit);
    Xbest  = X(bidx,:);
    Xworst = X(widx,:);

    for i = 1:N

        if fes >= MaxFEs
            best_f = min(fit);
            return;
        end

        r1 = rand(1,D);
        r2 = rand(1,D);

        % ---- standard JAYA update (NO abs) ----
        Xj = X(i,:) + r1.*(Xbest - X(i,:)) - r2.*(Xworst - X(i,:));
        Xj = min(max(Xj, lb), ub);

        fj = cec22_func(Xj, func_id);
        fes = fes + 1;

        if fj < fit(i)
            X(i,:) = Xj;
            fit(i) = fj;
        end
    end
end

best_f = min(fit);
end
