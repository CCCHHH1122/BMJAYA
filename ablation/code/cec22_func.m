function f = cec22_func(X, func_id)
% CEC2022 evaluation wrapper for cec22_test_func.mexw64
% Supports:
%   - Single solution: D×1 or 1×D
%   - Population: N×D  (each row is a solution)
%
% Returns:
%   - Scalar for single solution
%   - N×1 vector for population

if isvector(X)
    % single solution
    x = X(:);                 % force D×1
    f = cec22_test_func(x, func_id);
else
    % population N×D
    N = size(X,1);
    f = zeros(N,1);
    for i = 1:N
        f(i) = cec22_test_func(X(i,:).', func_id);  % D×1
    end
end
end
