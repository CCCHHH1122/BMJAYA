# BMJAYA

This repository contains only the standalone BMJAYA optimizer core.

BMJAYA is a parameter-free multi-strategy enhanced JAYA algorithm for
global optimization. The implementation keeps the core mechanisms used in
the manuscript:

- multipopulation cooperative evolution;
- JAYA update attempted first for each individual;
- conditional Best-Mean-Random repair when the JAYA update fails;
- global elite sharing among subpopulations;
- fixed internal subpopulation count `K = 4`.

## Files

- `BMJAYA.m`  
  Standalone MATLAB implementation of BMJAYA. It accepts a general
  objective function handle and does not depend on CEC benchmark files.

- `examples/demo_sphere.m`  
  Minimal example on the Sphere function.

## Usage

```matlab
D = 30;
N = 50;
MaxFEs = 100000;
lb = -100;
ub = 100;

objfun = @(x) sum(x.^2);

[best_f, best_x, curve] = BMJAYA(objfun, D, N, MaxFEs, lb, ub);
```
