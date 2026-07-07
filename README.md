# BMJAYA

This repository contains MATLAB source code, benchmark settings, and experimental result files for the manuscript:

**BMJAYA: A parameter-free multi-strategy enhanced JAYA algorithm for global optimization**

BMJAYA enhances the original JAYA framework by integrating a conditionally triggered Best--Mean--Random (BMR) repair mechanism, multipopulation cooperative evolution, and global elite sharing. The fixed internal subpopulation setting used in the reported experiments is `M = 4`. This setting is used as a structural organization choice and is not tuned separately for individual benchmark functions.

## Repository Structure

- `cec2022/code`
  MATLAB code for the CEC-2022 benchmark experiments with PSO, GWO, DE, WSO, JAYA, SBOA, DOA, and BMJAYA.

- `cec2022/results`
  Final CEC-2022 result CSV files used for the 10D and 20D comparison tables.

- `cec2022/figures`
  Representative convergence and radar-chart figures.

- `cec2017_50D/code`
  MATLAB code for the CEC-2017 50D validation experiments.

- `cec2017_50D/results`
  Final CEC-2017 50D result CSV files and convergence figure.

- `ablation/code`
  MATLAB code for the ablation study and subpopulation-sensitivity experiments.

- `ablation/results`
  Result CSV files for the ablation study and K-sensitivity analysis.

- `engineering_validation`
  Notes for the BMJAYA-XGBoost engineering validation. The corrosion dataset is available from the corresponding author upon reasonable request.

## Main Experimental Settings

- Population size: `N = 50`
- Independent runs: `30`
- CEC-2022 dimensions: `D = 10` and `D = 20`
- CEC-2022 MaxFEs: `100000`
- CEC-2017 dimension: `D = 50`
- CEC-2017 MaxFEs: `500000`
- BMJAYA fixed internal subpopulation setting: `M = 4`

## Re-running CEC-2022 Experiments

In MATLAB:

```matlab
cd('cec2022/code')
run_cec2022_all
```

The CEC-2022 support files and `input_data` folder must remain in `cec2022/code`.

## Re-running CEC-2017 50D Experiments

In MATLAB:

```matlab
cd('cec2017_50D/code')
run_cec2017_50D_all
```

The CEC-2017 support files and `input_data` folder must remain in `cec2017_50D/code`.

## Re-running the Ablation Study

In MATLAB:

```matlab
cd('ablation/code')
run_cec2022_all
```

For K-sensitivity:

```matlab
run_k_sensitivity_mjaya_bmjaya
```

## Notes

- The CEC benchmark MEX files are included for the Windows MATLAB environment used in the experiments.
- If another operating system is used, the corresponding CEC benchmark MEX files should be compiled from the official CEC source code.
- The repository intentionally does not include SFOA, because SFOA was not used as a comparison algorithm in the revised manuscript.
- The engineering corrosion dataset is not redistributed here and is available from the corresponding author upon reasonable request.

