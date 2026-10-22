# M3a — partially pooled group coefficients

[Model comparison](../README.md) · [Project overview](../../README.md)

I fitted M3a on the external public cohort and retained three CV fold fits, each with four chains. M3a is an alternative fitted specification, not the packaged M3b deployment mixture.

[Executable Stan definition](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/workflow/stan/m3a.stan) · [Architecture diagram](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m3a/figures/architecture.svg)

## Inputs and model boundary

`X` contains candidate features and `group_id` assigns each candidate row to a fitted group. The meaning and mapping of the grouping factor must come from the training metadata; the Stan file does not define a disease or mechanism taxonomy.

M3a uses a non-centered coefficient hierarchy. It does not add a copula contribution. The shared pseudo-observation declarations and helper copula function do not enter this predictor.

## Non-centered hierarchy

For group g and feature k:

```math
\beta_{gk}=\beta_k+\tau_{\beta,k}z_{\beta,gk}.
```

The coefficient vector for a group is centered on a shared population vector, with feature-specific scales controlling variation between groups.

```math
\beta_k\sim\mathcal{N}(0,s_{\beta,k}^2),\qquad
z_{\beta,gk}\sim\mathcal{N}(0,1).
```

```math
\tau_{\beta,k}\sim\mathrm{HalfNormal}(s_{\tau\beta,k}).
```

The scale inputs are `beta_prior_scale` and `tau_beta_prior_scale`. The Stan model constrains `tau_beta` to be nonnegative and includes the half-Normal normalization term K log(2). These constants are part of the density specification, especially when comparing marginal likelihoods.

## Candidate score and likelihood

```math
\eta_i=\sum_k x_{ik}\beta_{g(i),k}.
```

For case c with its supplied candidate set:

```math
p_{ci}(\theta)=\frac{\exp(\eta_i(\theta))}
{\sum_{j\in I_c}\exp(\eta_j(\theta))},\qquad
\log L(\theta)=\sum_c\log p_{c,y_c}(\theta).
```

`y_pos` identifies the labeled candidate locally within each contiguous case segment. Group-specific coefficients share population information; they are not unrelated regressions fitted independently for each group.

## Prediction for known and unseen groups

For an already fitted group, generated quantities use its posterior group coefficients. The `group_new = G + 1` slot represents an unseen group and receives a fresh population-derived vector at each posterior draw:

```math
\beta_{*,k}^{(s)}=\beta_k^{(s)}+
\tau_{\beta,k}^{(s)}\epsilon_k^{(s)},\qquad
\epsilon_k^{(s)}\sim\mathcal{N}(0,1).
```

All rows assigned to that unseen-group slot share the sampled vector within the draw. The implementation does not draw an independent group effect for every candidate row or represent multiple distinct unseen groups with separate slots.

Generated quantities return `p_new` after within-case softmax. Average probabilities across draws, not coefficients before softmax. These are candidate-set-conditional probabilities, not clinical pathogenicity probabilities.

## Retained-fit diagnostics

[Diagnostic table](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m3a/figures/m3a_diagnostics.png)

- [Fold 1: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m3a/folds/fold01/README.md)
- [Fold 2: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m3a/folds/fold02/README.md)
- [Fold 3: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m3a/folds/fold03/README.md)

I report retained-parameter diagnostics and recorded divergence counts separately. Treedepth remains unverified. The table percentage is recorded divergences divided by retained draws, not an independently verified transition-rate audit. Parameter diagnostics do not establish predictive performance.

## Implementation and artifacts

This page is on `main`; its linked implementation and diagnostic figures are currently on `external-training`. Training/scoring entrypoints and shared dependencies remain under integration. This page does not create or validate a complete runner.

[Retained-fold release](https://github.com/LukefromWA/rare-disease-evidence-pipeline/releases/tag/retained-folds-v1): availability depends on publication and repository access. Artifact identities and fingerprints are recorded on the fold pages. No retained model file is changed here.
