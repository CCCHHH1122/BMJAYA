# CEC-2022 Experiments

This folder contains the code and results for the CEC-2022 benchmark experiments.

Run from MATLAB:

```matlab
cd('cec2022/code')
run_cec2022_all
```

The compared algorithms are PSO, GWO, DE, WSO, JAYA, SBOA, DOA, and BMJAYA. All algorithms use the same population size, maximum number of function evaluations, and number of independent runs.

BMJAYA uses the fixed internal subpopulation setting `M = 4`.

