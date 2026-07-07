function f = cec22_func(X, func_id)
% CEC2022 evaluation wrapper for cec22_test_func.mexw64.
% Supports a single solution or a population stored row-wise.

global CEC22_CURVE_RECORDER

if isvector(X)
    x = X(:);
    f = cec22_test_func(x, func_id);
else
    N = size(X, 1);
    f = zeros(N, 1);
    for i = 1:N
        f(i) = cec22_test_func(X(i, :).', func_id);
    end
end

if isstruct(CEC22_CURVE_RECORDER) ...
        && isfield(CEC22_CURVE_RECORDER, 'active') ...
        && CEC22_CURVE_RECORDER.active
    values = f(:);
    n = numel(values);
    start_idx = CEC22_CURVE_RECORDER.fes + 1;
    end_idx = min(CEC22_CURVE_RECORDER.fes + n, CEC22_CURVE_RECORDER.MaxFEs);
    if start_idx <= end_idx
        keep_n = end_idx - start_idx + 1;
        CEC22_CURVE_RECORDER.values(start_idx:end_idx) = values(1:keep_n);
    end
    CEC22_CURVE_RECORDER.fes = CEC22_CURVE_RECORDER.fes + n;
end
end
