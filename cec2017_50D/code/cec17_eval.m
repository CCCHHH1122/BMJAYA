function f = cec17_eval(X, func_id)
% CEC17_EVAL Safe adapter for official cec17_func MEX.
%
% Input:
%   X       - one candidate (1xD or Dx1) OR population (N x D)
%   func_id - CEC2017 function ID
%
% Output:
%   f       - column vector of objective values
%
% Official cec17_func generally expects D x N, one candidate per column.

if isvector(X)
    f = cec17_func(X(:), func_id);
else
    f = cec17_func(X.', func_id);
end

f = f(:);

assert(all(isfinite(f)), ...
    'cec17_func returned NaN or Inf for F%d.', func_id);
end
