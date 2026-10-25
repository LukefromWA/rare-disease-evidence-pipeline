# Bayesian candidate-ranking models

[Project overview](../README.md) · [Workflow guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/workflow/README.md) · [Workflow library](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow/lib)

I fitted four Stan specifications on an external public cohort and retained three cross-validation fold fits per specification, with four chains per fit. The models compare global versus partially pooled feature coefficients, with or without an added Gaussian-copula contribution. The hackathon child's genome is not the training data for these fits.

| Model | Feature coefficients | Added dependence contribution | Recorded role |
|---|---|---|---|
| [M2a](m2a/README.md) | Global | None | Alternative fitted specification |
| [M2b](m2b/README.md) | Global | Global copula correlation and utility | Alternative fitted specification |
| [M3a](m3a/README.md) | Partially pooled by group | None | Alternative fitted specification |
| [M3b](m3b/README.md) | Partially pooled by group | Partially pooled correlation and utility | Three-component frozen deployment mixture |

## What these pages contain

Each page explains the inputs, score, likelihood, active priors, prediction behavior, interpretation and retained-fit qualifications. Equations describe the inspected Stan definitions, not a new model fitted by this documentation update. Display equations use fenced math blocks rather than bare LaTeX bracket delimiters, and avoid the unsupported operatorname macro.

These model descriptions are shared between `main` and `workflow`. The implementation, architecture figures and retained diagnostic material are currently linked from `external-training`; this documentation commit does not merge those files or create a runnable workflow or library. The workflow guide documents execution routes. Shared implementation files are maintained in `workflow/lib/` on the `workflow` branch; the top-level `lib/` directory is a separate consolidation placeholder.

## Shared interpretation

Each training case has one labeled candidate position in its supplied candidate set. All four specifications use within-case categorical-logit likelihoods and return within-case softmax probabilities for supplied new cases. Feature meanings, column order, group mapping and transformations come from the matching training record, not the Stan filename.

A probability conditional on a supplied candidate set is not a calibrated probability of clinical pathogenicity. Normalized A-layer scores, ranking-sensitivity frequencies and Stan probabilities must remain separately labeled.

M3b deployment averages three fold-specific posterior predictive distributions equally. It is not a single full-data posterior or a Bayesian model average across the four specifications. OOF evaluation must use the component that excluded the evaluated case.

## Diagnostics and reproduction

Parameter convergence summaries, sampler pathologies and held-out predictive evaluation answer different questions. Recorded divergences are reported separately; treedepth remains unverified. A percentage formed from recorded divergences divided by retained draws is not independently established as a same-attempt, unthinned post-warmup transition rate.

Training/scoring entrypoints and shared dependencies remain under integration. These pages do not claim an end-to-end reproduction from `main` or complete scientific/privacy clearance.

[Existing diagnostic overview](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/external-training/models) · [Stan definitions](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/external-training/workflow/stan) · [Retained-artifact documentation](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/external-training/docs/retained-artifacts.md)

Research use and human review only; not a diagnostic or treatment-selection system.
