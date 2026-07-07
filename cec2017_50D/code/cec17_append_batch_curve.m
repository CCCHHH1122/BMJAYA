function [curve, fes, best_f] = cec17_append_batch_curve(curve, fes, values, best_f)
% Append one or multiple FE values to an FE-based best-so-far curve.
values = values(:);
n = numel(values);

if n == 0
    return;
end

assert(fes + n <= numel(curve), ...
    'Attempted to exceed MaxFEs.');

batch_best = cummin(values);
batch_best = min(batch_best, best_f * ones(size(batch_best)));

curve(fes + (1:n)) = batch_best;
fes = fes + n;
best_f = batch_best(end);
end
