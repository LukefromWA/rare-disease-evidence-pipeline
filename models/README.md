# Retained models

Four Stan specifications were fitted, with three retained CV folds per specification. Model pages describe their statistical structure and diagnostics; fold pages show individual plots and identify the corresponding retained object.

- [M2a: global coefficients](m2a/README.md)
- [M2b: global coefficients and copula dependence](m2b/README.md)
- [M3a: partially pooled coefficients](m3a/README.md)
- [M3b: partially pooled coefficients and dependence](m3b/README.md)

The downstream retained predictor uses the three M3b fold components with equal weights. Retaining all four specifications does not imply they were blended into a single predictor.

[Retained-artifact guide](../docs/retained-artifacts.md) · [Executable model definitions](../workflow/stan/README.md)
