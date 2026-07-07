CEC2017 50D strict-MaxFEs package for BMJAYA experiments
==========================================================

1. Active algorithms
--------------------
PSO, GWO, WSO, BMJAYA, DE, JAYA, DOA, SBOA

2. Required official benchmark files
------------------------------------
Download the official CEC2017 single-objective bound-constrained MATLAB
benchmark implementation and place these items under your CEC2017 folder:

  cec17_func.mexw64
  input_data\

Then update this line near the top of run_cec2017_50D_all.m:

  CEC_DIR = 'E:\桌面\project\cec2017';

3. Run benchmark
----------------
In MATLAB:

  run_cec2017_50D_all

The default experiment is:

  D       = 50
  N       = 50
  MaxFEs  = 100000
  Runs    = 30
  FuncIDs = [1, 3:30]

F2 is excluded by default because the original CEC2017 bound-constrained
benchmark implementation is commonly treated as unstable for F2.

4. Outputs
----------
The script creates results_cec2017\ and writes:

  CEC2017_<ALGORITHM>_50D.csv
  CEC2017_ALL_Algorithms_50D_with_Rank.csv
  CEC2017_ALL_Algorithms_50D_with_Rank_scientific.csv
  CEC2017_FinalRank_50D.csv
  results_CEC2017_8alg_D50_N50_FEs100000_Runs30.mat

The summary table uses RAW objective values f(x).
The convergence curves use |f(x)-f*| + 1e-12.

5. Plot convergence curves
--------------------------
After running the benchmark:

  plot_cec2017_50D_convergence

It exports:

  CEC2017_convergence_8alg_D50_median.png
  CEC2017_convergence_8alg_D50_median.pdf

6. Important notes
------------------
- DOA is a paper-derived adapter based on the equations in the supplied PDF.
  Validate it against the authors' official source code before presenting it
  as an exact official implementation.
- SBOA is adapted to strict MaxFEs.
- BMJAYA preserves your original BMR rule:
    Xbest + (Xmean-Xi) + (Xrand-Xi)
  with self-selection allowed and elite sharing every generation.
