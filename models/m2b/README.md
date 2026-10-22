# M2b — retained external-training fits

I report the three retained CV fold fits for M2b here. Each fit has four chains.

[Stan model source](../../workflow/stan/m2b.stan)

## Fold diagnostic table

![M2b retained fold diagnostics](figures/m2b_diagnostics.png)

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

![M2B model architecture](figures/architecture.svg)

### Global coefficients with a Gaussian-copula term

A masked Gaussian-copula log-density term augments the linear score. The copula correlation is global, as is its contribution coefficient.

\[
\eta_i = x_i^{\mathsf T}\beta + b_{\mathrm{cop}}\,m_i\log c_{\rho}(u_{H,i},u_{N,i}).
\]

Here, `X` is the prepared feature matrix. For the copula specifications, the HPO/network pseudo-observations are separate inputs, and `m_i` is the availability mask. Feature meanings and column order come from the matching preprocessing record, not the Stan filename.

### Likelihood and prediction

Training uses a categorical-logit likelihood over the candidate set within each case, with `y_pos` identifying the single labeled candidate. Generated quantities apply a within-case softmax to new candidate scores and return `p_new`. These probabilities are conditional on the supplied candidate set, not clinical pathogenicity probabilities.

### Parameters and priors

`beta`, `rho_raw`, and `beta_cop` have zero-centered Normal priors. The correlation is `tanh(rho_raw)`. Prior scales are supplied as data.

[Read the executable Stan definition](../../workflow/stan/m2b.stan).

### Retained-fit downloads

[Open the retained-fold release](https://github.com/LukefromWA/rare-disease-evidence-pipeline/releases/tag/retained-folds-v1). Individual artifact links and fingerprints appear on the fold pages.
<!-- END MODEL_ARCHITECTURE -->
