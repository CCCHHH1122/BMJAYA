# Ablation Study

This folder contains the code and results used for the BMJAYA ablation study.

The ablation variants are:

- `JAYA`: original JAYA with strict MaxFEs
- `BJAYA`: JAYA with the BMR repair mechanism
- `MJAYA`: multipopulation JAYA with elite sharing
- `BMJAYA`: BMR repair + multipopulation evolution + elite sharing

Run from MATLAB:

```matlab
cd('ablation/code')
run_cec2022_all
```

For K-sensitivity:

```matlab
run_k_sensitivity_mjaya_bmjaya
```

The manuscript reports BMJAYA using the fixed internal subpopulation setting `M = 4`.

