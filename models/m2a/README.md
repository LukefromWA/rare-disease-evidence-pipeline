# M2a — global conditional-softmax model

[Model comparison](../README.md) · [Project overview](../../README.md)

I fitted M2a on the external public cohort and retained three CV fold fits, each with four chains. M2a is an alternative fitted specification, not the packaged M3b deployment mixture.

[Executable Stan definition](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/workflow/stan/m2a.stan) · [Architecture diagram](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m2a/figures/architecture.svg)

## Inputs and model boundary

`X` is the prepared candidate-by-feature matrix. `case_start` and `case_end` define contiguous candidate segments; `y_pos` is the local index of the single labeled candidate in each training case. Column meanings, feature transformations and prior scales must match the retained run record.

M2a uses one coefficient vector across all groups. Although the shared Stan data interface declares group labels and HPO/network pseudo-observations, they do not enter M2a's score. The helper copula function is also not used by this predictor. Omitting an added copula term does not imply that the input features are empirically independent.

## Global coefficients

For candidate feature row i, the model score is:

```math
\eta_i=x_i^{\mathsf T}\beta.
```

The active coefficient priors are:

```math
\beta_k\sim\mathcal{N}(0,s_{\beta,k}^2).
```

The scale vector is supplied as `beta_prior_scale`. The other shared prior-scale declarations are not additional active priors in M2a.

## Case-conditional likelihood

For the candidate rows of case c, with labeled candidate y_c:

```math
p_{ci}(\beta)=\frac{\exp(\eta_i)}{\sum_{j\in I_c}\exp(\eta_j)}.
```

```math
\log L(\beta)=\sum_c\left[\eta_{y_c}-\log\sum_{j\in I_c}\exp(\eta_j)\right].
```

In Stan, the observed label is passed as its local `y_pos` within the case segment, not a global matrix-row index. Case membership must therefore be preserved when constructing and packing the data.

## New-case prediction

Generated quantities use `X_new` with the sampled global coefficients and apply softmax separately within each new case, returning `p_new`. No refitting is needed to generate these conditional predictions from retained posterior draws.

```math
\widehat p_{ci}=\frac1S\sum_{s=1}^{S}
\frac{\exp(x_i^{\mathsf T}\beta^{(s)})}
{\sum_{j\in I_c}\exp(x_j^{\mathsf T}\beta^{(s)})}.
```

Average probabilities across draws after the within-case softmax; softmax of an averaged score is not the same quantity. Predictions depend on the supplied candidate set and are not clinical pathogenicity probabilities.

## Retained-fit diagnostics

[Diagnostic table](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m2a/figures/m2a_diagnostics.png)

- [Fold 1: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m2a/folds/fold01/README.md)
- [Fold 2: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m2a/folds/fold02/README.md)
- [Fold 3: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m2a/folds/fold03/README.md)

I report retained-parameter diagnostics and recorded divergence counts separately. Treedepth remains unverified. The table percentage is recorded divergences divided by retained draws, not an independently verified transition-rate audit. Parameter diagnostics do not establish predictive performance.

## Implementation and artifacts

The integrated workflow implementation is maintained on the `workflow` branch. The linked Stan definitions and retained diagnostic figures are maintained on `external-training`. Training/scoring entrypoints and their shared dependencies remain under integration. This page does not create those components or certify an end-to-end reproduction.

[Retained-fold release](https://github.com/LukefromWA/rare-disease-evidence-pipeline/releases/tag/retained-folds-v1): availability depends on release publication and repository access. The fold documentation records artifact identities and fingerprints. No retained model file is changed by this documentation update.
