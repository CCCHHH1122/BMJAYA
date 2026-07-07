function [lb, ub] = cec17_bounds(D)
% CEC17_BOUNDS Standard CEC2017 bound-constrained search domain.
lb = -100 * ones(1, D);
ub =  100 * ones(1, D);
end
