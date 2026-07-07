function [best_f, curve, fes] = alg_wso_compact_cec17(func_id, D, N, MaxFEs, lb, ub)
% Compact WSO baseline with strict MaxFEs.
pmin = 0.5; pmax = 1.5;
tau = 4.125;
a0 = 6.25; a1 = 100; a2 = 0.0005;
fmin = 0.07; fmax = 0.75;

disc = max(tau^2 - 4 * tau, 0);
mu = (2 - tau - sqrt(disc)) / 2;

curve = nan(MaxFEs, 1);
fes = 0;

W = rand(N, D) .* (ub - lb) + lb;
V = zeros(N, D);

fit = cec17_eval(W, func_id);
[curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fit, inf);

P = W;
Pfit = fit;
[~, gidx] = min(Pfit);
G = P(gidx, :);

max_steps = max(1, ceil((MaxFEs - fes) / N));
step = 0;

while fes < MaxFEs
    step = step + 1;

    decay = exp(-(4 * step / max_steps)^2);
    p1 = pmax + (pmax - pmin) * decay;
    p2 = pmin + (pmax - pmin) * decay;

    f = fmin + (fmax - fmin) * step / max_steps;
    mv = 1 / (a0 + exp((max_steps / 2 - step) / a1));
    ss = abs(1 - exp(-a2 * step / max_steps));

    peer_idx = randi(N, N, D);
    lin = sub2ind([N, D], peer_idx, repmat(1:D, N, 1));
    Wvbest = P(lin);

    c1 = rand(N, D);
    c2 = rand(N, D);

    Vcand = mu * (V + p1 * (G - W) .* c1 + p2 * (Wvbest - W) .* c2);
    Wnew = W + Vcand ./ max(f, 1e-12);

    jump = rand(N, 1) < mv;
    if any(jump)
        Wnew(jump, :) = rand(sum(jump), D) .* (ub - lb) + lb;
    end

    Wnew = min(max(Wnew, lb), ub);

    r1 = rand(N, 1);
    r2 = rand(N, 1);
    r3 = rand(N, 1);

    Dw = abs(rand(N, D) .* (G - Wnew));
    dir = ones(N, 1);
    dir(r2 < 0.5) = -1;

    What = Wnew;
    mask = r3 <= ss;

    if any(mask)
        What(mask, :) = G + (r1(mask) .* dir(mask)) .* Dw(mask, :);
    end

    denom = 2 * rand(N, D);
    denom(denom == 0) = 1e-12;

    Wcand = min(max((Wnew + What) ./ denom, lb), ub);

    nEval = min(N, MaxFEs - fes);
    fnew = cec17_eval(Wcand(1:nEval, :), func_id);
    [curve, fes, best_f] = cec17_append_batch_curve(curve, fes, fnew, best_f);

    W(1:nEval, :) = Wcand(1:nEval, :);
    V(1:nEval, :) = Vcand(1:nEval, :);

    improved = fnew < Pfit(1:nEval);
    idx = find(improved);

    if ~isempty(idx)
        P(idx, :) = W(idx, :);
        Pfit(idx) = fnew(improved);
    end

    [pbest, gidx] = min(Pfit);
    G = P(gidx, :);
    best_f = min(best_f, pbest);
end

curve = cec17_fill_tail(curve);
end
