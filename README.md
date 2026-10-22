# External training — Bayesian model comparison

I trained four Stan specifications using an external public cohort and retained three cross-validation fits per specification. I use the three retained M3b components as an equal-weight frozen predictor for downstream inference.

## Retained-model diagnostics

![Retained model diagnostic table](figures/all_models_diagnostics.png)

I recomputed parameter diagnostics from the retained draws using their stored chain metadata. All twelve model-fold parameter checks passed the declared thresholds. I report the nonzero divergence counts separately, and treedepth remains unverified.

## Model documentation

| Model | Fold diagnostics and plots | Analysis role |
|---|---|---|
| M2a | [Model README](models/m2a/README.md) | Alternative fitted specification |
| M2b | [Model README](models/m2b/README.md) | Alternative fitted specification |
| M3a | [Model README](models/m3a/README.md) | Alternative fitted specification |
| M3b | [Model README](models/m3b/README.md) | Frozen downstream specification |

## Reading the diagnostics

I use R-hat and effective sample size to assess retained parameter draws. I do not treat these checks as predictive-performance comparisons. The divergence percentages use recorded counts divided by retained draws; the same-attempt and unthinned post-warmup assumptions must hold before calling them verified transition rates.

## Evaluation and frozen inference

I distinguish OOF evaluation from new-case scoring. An OOF case is evaluated with the component that excluded it from training. My frozen M3b predictor combines three components for new cases; it is not the identical OOF predictor or one full-data posterior. Each component has four chains.

## Reproduction status

I have included model definitions, diagnostic records and individual fold plots. I am still integrating and verifying the training/scoring entrypoints and their shared dependencies, so this checkout is not yet a complete end-to-end runner. I keep raw training objects and private-case predictions out of this source-and-documentation import.

## AI assistance

I used Perplexity Pro primarily for documentation, debugging and code/workflow integration. I configured the relevant privacy and data-usage settings from the beginning to prevent private data from being used for model training, and I checked those settings throughout the project. I enabled equivalent training-use restrictions from the outset on my other AI accounts.

I will identify any additional contributing providers, products and tiers in my final methods record. This describes the settings I verified; I am not claiming zero retention.

See my [AI assistance disclosure](docs/AI_ASSISTANCE.md).
