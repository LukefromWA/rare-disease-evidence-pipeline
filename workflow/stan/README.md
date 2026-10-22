# Stan model definitions

| Source | Coefficients | Dependence term | Documentation |
|---|---|---|---|
| [m2a.stan](m2a.stan) | Global | None in the score | [M2a](../../models/m2a/README.md) |
| [m2b.stan](m2b.stan) | Global | Global Gaussian copula | [M2b](../../models/m2b/README.md) |
| [m3a.stan](m3a.stan) | Partially pooled by group | None in the score | [M3a](../../models/m3a/README.md) |
| [m3b.stan](m3b.stan) | Partially pooled by group | Group-specific Gaussian copula | [M3b](../../models/m3b/README.md) |

All specifications use a categorical-logit likelihood within each training case and generate new-case candidate probabilities through a within-case softmax. Shared data declarations and helper functions do not imply every model uses every input in its score. Prior scales are supplied as data.
