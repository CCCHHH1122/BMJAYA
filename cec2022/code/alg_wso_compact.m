function best_f = alg_wso_compact(func_id, D, N, iters, lb, ub)

% Common WSO parameter defaults (baseline-like)
pmin = 0.5; pmax = 1.5;
tau  = 4.125;
a0   = 6.25; a1 = 100; a2 = 0.0005;
fmin = 0.07; fmax = 0.75;

disc = tau^2 - 4*tau;
disc = max(disc, 0);
mu = (2 - tau - sqrt(disc))/2;

W = rand(N,D).*(ub-lb) + lb;  % positions
V = zeros(N,D);               % velocities

fit = cec22_func(W, func_id);
P = W; Pfit = fit;

[best_f, gidx] = min(fit);
G = W(gidx,:);

K = max(iters, 1);

for k = 1:iters
    decay = exp(- (4*k/K)^2 );
    p1 = pmax + (pmax - pmin)*decay;
    p2 = pmin + (pmax - pmin)*decay;

    f = fmin + (fmax - fmin)*(k/K);

    mv = (1/(a0 + exp((K/2 - k)/a1)));  % monotone increasing-ish
    ss = abs(1 - exp(-a2*k/K));

    % select peer-best positions (vectorized sampling)
    peer_idx = randi(N, N, D);
    lin = sub2ind([N,D], peer_idx, repmat(1:D, N, 1));
    Wvbest = P(lin);

    c1 = rand(N,D);
    c2 = rand(N,D);

    V = mu*(V + p1*(G - W).*c1 + p2*(Wvbest - W).*c2);

    % move toward prey
    Wnew = W + V./max(f,1e-12);

    % exploration jump with probability mv (simple stable variant)
    jump = rand(N,1) < mv;
    if any(jump)
        Wnew(jump,:) = rand(sum(jump),D).*(ub-lb) + lb;
    end

    Wnew = min(max(Wnew, lb), ub);

    % move toward best shark (simplified eq.12)
    r1 = rand(N,1); r2 = rand(N,1); r3 = rand(N,1);
    Dw = abs(rand(N,D).*(G - Wnew));
    dir = ones(N,1); dir(r2 < 0.5) = -1;

    What = Wnew;
    mask = (r3 <= ss);
    if any(mask)
        What(mask,:) = G + (r1(mask).*dir(mask)).*Dw(mask,:);
    end

    denom = 2*rand(N,D);
    denom(denom==0) = 1e-12;
    W = (Wnew + What)./denom;
    W = min(max(W, lb), ub);

    fit = cec22_func(W, func_id);

    improved = fit < Pfit;
    P(improved,:) = W(improved,:);
    Pfit(improved) = fit(improved);

    [cur_best, gidx] = min(fit);
    if cur_best < best_f
        best_f = cur_best;
        G = W(gidx,:);
    end
end

end
