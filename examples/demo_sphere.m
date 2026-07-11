clear; clc;

D = 30;
N = 50;
MaxFEs = 100000;
lb = -100;
ub = 100;

objfun = @(x) sum(x.^2);

[best_f, best_x, curve] = BMJAYA(objfun, D, N, MaxFEs, lb, ub);

fprintf('Best objective value: %.12e\n', best_f);
fprintf('Best solution norm:   %.12e\n', norm(best_x));

figure;
semilogy(curve, 'LineWidth', 1.5);
grid on;
xlabel('Function evaluations');
ylabel('Best-so-far objective value');
title('BMJAYA convergence on Sphere');
