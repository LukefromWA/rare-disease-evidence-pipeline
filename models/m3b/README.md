# M3b — partially pooled coefficients and copula dependence

[Model comparison](../README.md) · [Project overview](../../README.md)

I fitted M3b on the external public cohort and retained three CV fold fits, each with four chains. The frozen downstream predictor combines those three M3b components with equal weights. The hackathon child's genome is not the training data for these fits.

[Executable Stan definition](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/workflow/stan/m3b.stan) · [Architecture diagram](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m3b/figures/architecture.svg)

## Inputs and model boundary

`X` supplies prepared candidate features; `group_id` supplies their group mapping. HPO/network pseudo-observations are separate inputs. The `copula_available` mask controls the added dependence contribution. Feature meanings, group definitions and transformations must match the retained run record.

M3b partially pools feature coefficients, copula correlation and copula utility. It does not impose a single global dependence contribution on every group, and the fitted copula term does not establish causal dependence between evidence sources.

## Feature-coefficient hierarchy

```math
\beta_{gk}=\beta_k+\tau_{\beta,k}z_{\beta,gk},\qquad
\beta_k\sim\mathcal{N}(0,s_{\beta,k}^2).
```

```math
z_{\beta,gk}\sim\mathcal{N}(0,1),\qquad
\tau_{\beta,k}\sim\mathrm{HalfNormal}(s_{\tau\beta,k}).
```

The coefficient scale inputs are `beta_prior_scale` and `tau_beta_prior_scale`.

## Dependence and copula-utility hierarchies

In the equations, the population correlation location corresponds to `mu_rho_raw`; the copula utility location, scale and standardized offsets correspond to `beta_cop`, `tau_cop` and `z_cop`.

```math
\rho_g=\tanh(\mu_\rho+\tau_\rho z_{\rho,g}),\qquad
\gamma_g=\gamma+\tau_\gamma z_{\gamma,g}.
```

The latent group offsets have independent standard-Normal priors conditional on the population parameters:

```math
z_{\rho,g}\sim\mathcal{N}(0,1),\qquad
z_{\gamma,g}\sim\mathcal{N}(0,1).
```

```math
\mu_\rho\sim\mathcal{N}(0,s_{\mu\rho}^2),\qquad
\gamma\sim\mathcal{N}(0,s_\gamma^2).
```

```math
\tau_\rho\sim\mathrm{HalfNormal}(s_{\tau\rho}),\qquad
\tau_\gamma\sim\mathrm{HalfNormal}(s_{\tau\gamma}).
```

The scale inputs are `mu_rho_prior_scale`, `tau_rho_prior_scale`, `beta_copula_prior_scale`, and `tau_copula_prior_scale`. The model includes the half-Normal normalization constants. The copula utility coefficients may be negative.

## Gaussian-copula contribution and score

```math
z_{hi}=\Phi^{-1}(u_{hi}),\qquad z_{ni}=\Phi^{-1}(u_{ni}).
```

```math
\ell_i(\rho)=-\frac12\log(1-\rho^2)
-\frac{\rho^2(z_{hi}^2+z_{ni}^2)-2\rho z_{hi}z_{ni}}{2(1-\rho^2)}.
```

```math
\eta_i=\sum_k x_{ik}\beta_{g(i),k}
+\gamma_{g(i)}m_i\ell_i(\rho_{g(i)}).
```

The availability mask removes the copula contribution when unavailable; it does not remove the linear features. Pseudo-observations must still satisfy the Stan input bounds because their inverse-Normal transforms are computed before masking.

## Case-conditional likelihood

For the supplied candidate rows of case c:

```math
p_{ci}(\theta)=\frac{\exp(\eta_i(\theta))}
{\sum_{j\in I_c}\exp(\eta_j(\theta))},\qquad
\log L(\theta)=\sum_c\log p_{c,y_c}(\theta).
```

`y_pos` is the local position of the single labeled candidate in its contiguous case segment. Generated quantities return `p_new` after applying the corresponding score and softmax within each new case. These probabilities depend on the supplied candidate set and are not clinical pathogenicity probabilities.

## New-group prediction

At each posterior draw, the `G + 1` slot samples a new coefficient vector, correlation and copula utility from the fitted population distributions:

```math
\beta_{*,k}=\beta_k+\tau_{\beta,k}\epsilon_{\beta,k},\qquad
\rho_*=\tanh(\mu_\rho+\tau_\rho\epsilon_\rho),\qquad
\gamma_*=\gamma+\tau_\gamma\epsilon_\gamma.
```

The epsilon terms are standard-Normal draws. Rows assigned to that unseen-group slot share its sampled parameters within a draw; the implementation does not allocate independent unseen-group effects to every row.

## Frozen three-component predictor and OOF distinction

For a new case, average the three component posterior predictive probabilities equally:

```math
\widehat p_{ci}^{\mathrm{frozen}}=
\frac13\sum_{f=1}^{3}\left[
\frac1{S_f}\sum_{s=1}^{S_f}p_{ci}(\theta_f^{(s)})
\right].
```

Each component contributes one third regardless of its retained draw count. Apply within-case softmax before averaging probabilities across draws and components.

This is a mixture of fold-specific posterior predictions, not a single full-data posterior and not evidence-weighted Bayesian model averaging across M2a, M2b, M3a and M3b. OOF evaluation uses the component that excluded the evaluated case. Applying the new-case mixture to a training-cohort case is not the same OOF procedure.

## Retained-fit diagnostics

[Diagnostic table](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m3b/figures/m3b_diagnostics.png)

- [Fold 1: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m3b/folds/fold01/README.md)
- [Fold 2: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m3b/folds/fold02/README.md)
- [Fold 3: trace, rank and interval plots](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/models/m3b/folds/fold03/README.md)

I report retained-parameter diagnostics and recorded divergence counts separately. Treedepth remains unverified. The table percentage is recorded divergences divided by retained draws, not an independently verified transition-rate audit. Parameter diagnostics do not establish predictive performance.

## Implementation and artifacts

This model description is included on the `workflow` branch; its linked Stan implementation and diagnostic figures are maintained on `external-training`. The documented deployment role does not imply the complete tested scoring adapter and its dependencies are already imported into `main`.

[Retained-fold release](https://github.com/LukefromWA/rare-disease-evidence-pipeline/releases/tag/retained-folds-v1): availability depends on publication and repository access. The fold pages record artifact identities and fingerprints. This documentation update does not change any fit, predictive weight, Stan definition, release asset or repository visibility setting.
