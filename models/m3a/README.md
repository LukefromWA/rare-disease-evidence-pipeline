# M3a — retained external-training fits

I report the three retained CV fold fits for M3a here. Each fit has four chains.

[Stan model source](../../workflow/stan/m3a.stan)

## Fold diagnostic table

![M3a retained fold diagnostics](figures/m3a_diagnostics.png)

I report the retained-parameter diagnostics and recorded divergence counts separately. Treedepth remains unverified. The table percentage is recorded divergences divided by retained draws, not an independently verified transition-rate audit.

## Individual fold plots

- [Fold 1: trace, rank and interval plots](folds/fold01/README.md)
- [Fold 2: trace, rank and interval plots](folds/fold02/README.md)
- [Fold 3: trace, rank and interval plots](folds/fold03/README.md)

## Analysis role

I fitted this alternative specification; I did not package it as the M3b deployment mixture.

## Reproduction status

I have included the Stan definition and diagnostic outputs. I am still integrating the tested training/scoring entrypoints and their shared dependencies into this checkout.

<!-- BEGIN MODEL_ARCHITECTURE -->
## Model structure

![M3A model architecture](figures/architecture.svg)

### Partially pooled group coefficients

Group-specific coefficients share population parameters through a non-centered hierarchy. There is no copula contribution in the score.

\[
\beta_g = \beta + \tau_{\beta}\odot z_{\beta,g}, \qquad \eta_i = x_i^{\mathsf T}\beta_{g(i)}.
\]

Here, `X` is the prepared feature matrix. For the copula specifications, the HPO/network pseudo-observations are separate inputs, and `m_i` is the availability mask. Feature meanings and column order come from the matching preprocessing record, not the Stan filename.

### Likelihood and prediction

Training uses a categorical-logit likelihood over the candidate set within each case, with `y_pos` identifying the single labeled candidate. Generated quantities apply a within-case softmax to new candidate scores and return `p_new`. These probabilities are conditional on the supplied candidate set, not clinical pathogenicity probabilities.

### Parameters and priors

`beta` has a zero-centered Normal prior; `z_beta` has standard Normal priors; nonnegative `tau_beta` has half-Normal priors. The implementation includes the half-Normal normalization constants.

[Read the executable Stan definition](../../workflow/stan/m3a.stan).

For a new group, generated quantities draw new group effects from the population distribution. The retained run metadata determines which groups were actually fitted.

### Retained-fit downloads

[Open the retained-fold release](https://github.com/LukefromWA/rare-disease-evidence-pipeline/releases/tag/retained-folds-v1). Individual artifact links and fingerprints appear on the fold pages.
<!-- END MODEL_ARCHITECTURE -->
