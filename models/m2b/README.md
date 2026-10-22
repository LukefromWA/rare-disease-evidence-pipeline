# M2b — global coefficients with a Gaussian-copula contribution

[Model comparison](../README.md) · [Project overview](../../README.md)

I fitted M2b on the external public cohort and retained three CV fold fits, each with four chains. M2b is an alternative fitted specification, not the packaged M3b deployment mixture.

[Executable Stan definition](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/workflow/stan/m2b.stan) · [Architecture diagram](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m2b/figures/architecture.svg)

## Inputs and model boundary

`X` supplies the prepared candidate features. HPO/network pseudo-observations `u_hpo` and `u_network` are separate inputs. The binary `copula_available` mask controls whether their copula contribution enters each row's score. All feature meanings and transformations must match the retained preprocessing record.

The coefficients, dependence correlation and copula utility coefficient are global, not group-specific. A copula term is a predictive dependence contribution; it is not evidence of a causal relationship between evidence sources.

## Gaussian-copula log density

Transform the supplied pseudo-observations:

```math
z_{hi}=\Phi^{-1}(u_{hi}),\qquad z_{ni}=\Phi^{-1}(u_{ni}),\qquad
\rho=\tanh(r).
```

Here r is `rho_raw`. The implemented log density is:

```math
\ell_i(\rho)=-\frac12\log(1-\rho^2)
-\frac{\rho^2(z_{hi}^2+z_{ni}^2)-2\rho z_{hi}z_{ni}}{2(1-\rho^2)}.
```

Stan bounds the pseudo-observations away from zero and one. They must remain valid finite inputs even for a masked row, because inverse-Normal transforms are computed before the contribution is masked.

## Candidate score and priors

Using gamma for `beta_cop` and m for the row's availability indicator:

```math
\eta_i=x_i^{\mathsf T}\beta+\gamma m_i\ell_i(\rho).
```

```math
\beta_k\sim\mathcal{N}(0,s_{\beta,k}^2),\qquad
r\sim\mathcal{N}(0,s_{\mu\rho}^2),\qquad
\gamma\sim\mathcal{N}(0,s_\gamma^2).
```

| Prior scale | Stan data field |
|---|---|
| Feature coefficients | `beta_prior_scale` |
| Raw correlation | `mu_rho_prior_scale` |
| Copula utility | `beta_copula_prior_scale` |

Gamma is not constrained to be nonnegative. A zero availability mask removes the added copula contribution, not the linear feature score. Shared hierarchy-scale data declarations do not create hierarchical parameters in M2b.

## Likelihood and prediction

For the candidate rows of case c:

```math
p_{ci}(\theta)=\frac{\exp(\eta_i(\theta))}
{\sum_{j\in I_c}\exp(\eta_j(\theta))},\qquad
\log L(\theta)=\sum_c\log p_{c,y_c}(\theta).
```

Stan's `y_pos` is the local index of the labeled candidate within its contiguous case segment. Generated quantities apply the same global score and within-case softmax to supplied new candidates, returning `p_new`.

Posterior predictive summaries average these probabilities across draws after softmax. The probabilities are conditional on the supplied candidate set and are not calibrated clinical pathogenicity probabilities.

## Retained-fit diagnostics

[Diagnostic table](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m2b/figures/m2b_diagnostics.png)

- [Fold 1: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m2b/folds/fold01/README.md)
- [Fold 2: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m2b/folds/fold02/README.md)
- [Fold 3: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m2b/folds/fold03/README.md)

I report retained-parameter diagnostics and recorded divergence counts separately. Treedepth remains unverified. The table percentage is recorded divergences divided by retained draws, not an independently verified transition-rate audit. Parameter diagnostics do not establish predictive performance.

## Implementation and artifacts

This page is on `main`; its linked implementation and diagnostic figures are currently on `external-training`. Training/scoring entrypoints and shared dependencies remain under integration. This documentation does not create or validate a full runner.

[Retained-fold release](https://github.com/LukefromWA/rare-disease-evidence-pipeline/releases/tag/retained-folds-v1): availability depends on publication and repository access. The fold pages record artifact identities and fingerprints. No retained model file is changed here.
