function best_f = alg_doa_wrapper(func_id, D, N, iters, lb, ub)
% ALG_DOA_WRAPPER
% Strict-budget Dhole Optimization Algorithm (DOA) adapter for the
% existing run_cec2022_all.m framework.
%
% Interface:
%   best_f = alg_doa_wrapper(func_id, D, N, iters, lb, ub)
%
% Runner convention:
%   iters  = floor((MaxFEs - N) / N)
%   MaxFEs = N + iters * N
%
% Output:
%   best_f = best RAW objective value f(x), not |f(x)-f*|.
%
% Objective function:
%   cec22_test_func(x(:), func_id)
%
% Reference:
%   Mohammed, B. O., Aghdasi, H. S., & Salehpour, P. (2025).
%   Dhole optimization algorithm: a new metaheuristic algorithm for
%   solving optimization problems. Cluster Computing, 28, 430.

assert(N >= 2, 'DOA requires N >= 2.');
assert(D >= 1, 'D must be positive.');
assert(iters >= 0, 'iters must be non-negative.');

lb = reshape(lb, 1, []);
ub = reshape(ub, 1, []);
if numel(lb) == 1, lb = repmat(lb, 1, D); end
if numel(ub) == 1, ub = repmat(ub, 1, D); end
assert(numel(lb) == D && numel(ub) == D, ...
    'lb and ub must be scalars or vectors with D elements.');
assert(all(lb < ub), ...
    'Every lower bound must be strictly smaller than the upper bound.');

MaxFEs = max(N, N + iters * N);
best_f = alg_doa_fixed_budget(func_id, D, N, MaxFEs, lb, ub);
end

% =====================================================================
function best_f = alg_doa_fixed_budget(func_id, D, N, MaxFEs, lb, ub)
% Strict-MaxFEs DOA implementation adapted from the paper equations.

% Paper parameter table
C1 = 1.0;
l  = 25.0;
k  = 0.5;
C3 = 3.0;
EF = 0.2;

% Initialization
X = rand(N, D) .* (ub - lb) + lb;
fit = inf(N, 1);
fes = 0;
for i = 1:N
    fit(i) = evaluate_cec22(X(i, :), func_id);
    fes = fes + 1;
end

[global_best_f, global_best_idx] = min(fit);
global_best_x = X(global_best_idx, :);

% One full generation evaluates N candidates.
T = max(1, ceil((MaxFEs - N) / N));
t = 0;

while fes < MaxFEs
    t = t + 1;

    [local_best_f, local_best_idx] = min(fit);
    local_best_x = X(local_best_idx, :);

    % Eq. (5)
    prey = (local_best_x + global_best_x) ./ 2;

    % Eq. (3): nominally 5..20 pack members
    PMN = round(rand * 15 + 5);

    % Eq. (7)
    C2 = 1 - min(t / T, 1);

    Xnew = X;

    for i = 1:N
        if fes >= MaxFEs
            best_f = global_best_f;
            return;
        end

        vocalization = rand;

        if vocalization < 0.5
            if PMN < 10
                % Eq. (6): searching stage
                Xcand = X(i, :) ...
                      + C2 .* rand(1, D) .* (prey - X(i, :));
            else
                % Eq. (8): encircling stage
                z = randi(N);
                while z == i
                    z = randi(N);
                end
                Xcand = X(i, :) - X(z, :) + prey;
            end
        else
            % Eq. (4): suitable time for hunting
            ps = C1 ./ (1 + exp(-k .* (PMN - l) ./ 2)) .* EF;

            % Eq. (10): prey size
            denominator = max(abs(local_best_f), eps);
            S = C3 .* rand .* abs(fit(i) ./ denominator);

            if S > (C3 + 1) ./ 2
                % Eq. (11): weakened prey
                Wprey = exp(-1 ./ max(S, eps)) .* local_best_x;

                % Eq. (12): repeated attack
                rcos = rand(1, D);
                rsin = rand(1, D);
                Xcand = X(i, :) ...
                      + Wprey .* ps .* ...
                        (cos(2 .* pi .* rcos) ...
                      - sin(2 .* pi .* rsin) .* Wprey .* ps);
            else
                % Eq. (13): immediate kill
                Xcand = (X(i, :) - global_best_x) .* ps ...
                      + ps .* rand(1, D) .* X(i, :);
            end
        end

        Xcand = min(max(Xcand, lb), ub);
        fcand = evaluate_cec22(Xcand, func_id);
        fes = fes + 1;

        % Paper pseudocode updates positions and fitness after each move.
        Xnew(i, :) = Xcand;
        fit(i) = fcand;

        if fcand < global_best_f
            global_best_f = fcand;
            global_best_x = Xcand;
        end
    end

    X = Xnew;
end

best_f = global_best_f;
end

% =====================================================================
function f = evaluate_cec22(x, func_id)
% Official CEC2022 MEX interface expects a D-by-1 column vector.
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
