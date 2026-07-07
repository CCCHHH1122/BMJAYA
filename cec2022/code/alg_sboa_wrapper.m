function best_f = alg_sboa_wrapper(func_id, D, N, iters, lb, ub)
% ALG_SBOA_WRAPPER
% Secretary Bird Optimization Algorithm (SBOA) wrapper for the existing
% iter-based CEC2022 experiment framework.
%
% Required calling interface:
%   best_f = algo_fn(func_id, D, N, iters, lb, ub)
%
% Compatible with:
%   iters = floor((MaxFEs - N) / N)
%
% Inverse mapping:
%   MaxFEs = N + iters * N
%
% This implementation:
%   1) Uses the official CEC2022 MEX function:
%         cec22_test_func(x(:), func_id)
%   2) Counts initialization evaluations.
%   3) Enforces the restored MaxFEs budget strictly.
%   4) Stops immediately when the FE budget is exhausted.
%   5) Returns the best RAW objective value f(x), not |f(x)-f*|.
%
% Put this file in the MATLAB project folder and add:
%   struct('name','SBOA','fn',@alg_sboa_wrapper)
% to the Algorithms cell array in run_cec2022_all.m.
%
% Reference:
%   Fu, Y., Liu, D., Chen, J., He, L. (2024).
%   Secretary bird optimization algorithm: a new metaheuristic for
%   solving global optimization problems.
%   Artificial Intelligence Review, 57, 123.
%
% Notes:
%   - This is a clean strict-budget MATLAB implementation adapted to your
%     current runner interface.
%   - The original SBOA contains hunting and predator-evasion phases.
%   - A complete SBOA cycle normally evaluates approximately 2*N
%     candidates after initialization.

% ------------------------- validation -------------------------
assert(N >= 3, 'SBOA requires N >= 3.');
assert(D >= 1, 'D must be positive.');
assert(iters >= 0, 'iters must be non-negative.');

lb = reshape(lb, 1, []);
ub = reshape(ub, 1, []);

if numel(lb) == 1
    lb = repmat(lb, 1, D);
end

if numel(ub) == 1
    ub = repmat(ub, 1, D);
end

assert(numel(lb) == D && numel(ub) == D, ...
    'lb and ub must be scalars or vectors with D elements.');

assert(all(lb < ub), ...
    'Every lower bound must be strictly smaller than the upper bound.');

% Restore the strict FE budget used by your existing experiment runner.
MaxFEs = max(N, N + iters * N);

best_f = alg_sboa_fixedN_budget(func_id, D, N, MaxFEs, lb, ub);
end

% =====================================================================
function best_f = alg_sboa_fixedN_budget(func_id, D, N, MaxFEs, lb, ub)
% ALG_SBOA_FIXEDN_BUDGET
% SBOA under a strict function-evaluation budget.
%
% Output:
%   best_f = best RAW objective value f(x)

% ------------------------- initialization -----------------------------
X = rand(N, D) .* (ub - lb) + lb;
fit = inf(N, 1);

fes = 0;
best_f = inf;
best_x = zeros(1, D);

for i = 1:N
    fit(i) = evaluate_cec22(X(i, :), func_id);
    fes = fes + 1;

    if fit(i) < best_f
        best_f = fit(i);
        best_x = X(i, :);
    end
end

% One complete SBOA cycle normally uses:
%   N evaluations in hunting phase
%   N evaluations in predator-evasion phase
%
% The stage schedule is therefore based on the number of complete
% FE-equivalent cycles allowed by the restored strict budget.
T_equiv = max(1, ceil((MaxFEs - N) / (2 * N)));
t = 0;

% ============================ main loop ================================
while fes < MaxFEs
    t = t + 1;
    ratio = min(t / T_equiv, 1);

    % ================================================================
    % Phase P1: hunting behavior
    % ================================================================
    for i = 1:N
        if fes >= MaxFEs
            best_f = min(best_f, min(fit));
            return;
        end

        if t < T_equiv / 3
            % Stage 1: searching for prey.
            idx = randperm(N, 2);
            Xcand = X(i, :) ...
                  + (X(idx(1), :) - X(idx(2), :)) .* rand(1, D);

        elseif t < 2 * T_equiv / 3
            % Stage 2: consuming prey using Brownian motion.
            RB = randn(1, D);
            Xcand = best_x ...
                  + exp(ratio^4) .* (RB - 0.5) .* (best_x - X(i, :));

        else
            % Stage 3: attacking prey using Levy flight.
            RL = 0.5 .* levy_flight(D);
            Xcand = best_x ...
                  + ((1 - ratio)^(2 * ratio)) .* X(i, :) .* RL;
        end

        Xcand = min(max(Xcand, lb), ub);
        fcand = evaluate_cec22(Xcand, func_id);
        fes = fes + 1;

        if fcand < fit(i)
            X(i, :) = Xcand;
            fit(i) = fcand;

            if fcand < best_f
                best_f = fcand;
                best_x = Xcand;
            end
        end
    end

    % ================================================================
    % Phase P2: escaping from predators
    % ================================================================
    for i = 1:N
        if fes >= MaxFEs
            best_f = min(best_f, min(fit));
            return;
        end

        if rand < 0.5
            % Strategy C1: camouflage by the environment.
            RB = randn(1, D);
            Xcand = best_x ...
                  + (2 .* RB - 1) .* ((1 - ratio)^2) .* X(i, :);

        else
            % Strategy C2: fly or run away.
            K = round(1 + rand);
            R2 = randn(1, D);
            random_idx = randi(N);
            Xrandom = X(random_idx, :);

            Xcand = X(i, :) ...
                  + R2 .* (Xrandom - K .* X(i, :));
        end

        Xcand = min(max(Xcand, lb), ub);
        fcand = evaluate_cec22(Xcand, func_id);
        fes = fes + 1;

        if fcand < fit(i)
            X(i, :) = Xcand;
            fit(i) = fcand;

            if fcand < best_f
                best_f = fcand;
                best_x = Xcand;
            end
        end
    end
end

best_f = min(best_f, min(fit));
end

% =====================================================================
function f = evaluate_cec22(x, func_id)
% Evaluate one SBOA candidate using the official CEC2022 MEX interface.
f = cec22_test_func(x(:), func_id);

assert(isscalar(f) && isfinite(f), ...
    'cec22_test_func returned a non-finite or non-scalar value.');

global CEC22_CURVE_RECORDER
if isstruct(CEC22_CURVE_RECORDER) ...
        && isfield(CEC22_CURVE_RECORDER, 'active') ...
        && CEC22_CURVE_RECORDER.active
    idx = CEC22_CURVE_RECORDER.fes + 1;
    if idx <= CEC22_CURVE_RECORDER.MaxFEs
        CEC22_CURVE_RECORDER.values(idx) = f;
    end
    CEC22_CURVE_RECORDER.fes = CEC22_CURVE_RECORDER.fes + 1;
end
end

% =====================================================================
function step = levy_flight(D)
% Mantegna Levy-flight step used by SBOA.
beta = 1.5;
scale = 0.01;

sigma = (gamma(1 + beta) * sin(pi * beta / 2) / ...
        (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);

u = randn(1, D) .* sigma;
v = randn(1, D);

step = scale .* u ./ (abs(v).^(1 / beta) + eps);
end
