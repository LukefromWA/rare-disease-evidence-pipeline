# M2a — retained external-training fits

[Stan model source](../../workflow/stan/m2a.stan)

## Diagnostic values

| Fold | Attempt | Chains | Draws | Max R-hat | Min bulk ESS | Min tail ESS | Parameter gate | Divergences | Count / draws |
|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| 1 | 2 | 4 | 6,000 | 1.001624 | 2721.1 | 2695.5 | PASS | 0 | 0.00% |
| 2 | 1 | 4 | 4,000 | 1.003650 | 1788.0 | 1986.8 | PASS | 0 | 0.00% |
| 3 | 3 | 4 | 8,000 | 1.001836 | 3529.7 | 3921.7 | PASS | 0 | 0.00% |

Parameter diagnostics are not predictive validation. Recorded divergences remain qualified, treedepth is unverified, and conflicting historical summaries are preserved separately.

Count / draws means recorded divergences divided by retained parameter draws. A verified post-warmup transition rate requires the same retained attempt and unthinned draws.

## Overview

![M2a diagnostics](figures/m2a_overview_panel.png)

## Fold 1

![M2a fold 1 panel](figures/m2a_fold01_panel.png)

[Trace](figures/m2a_fold01_trace.png) · [Rank](figures/m2a_fold01_rank.png) · [Intervals](figures/m2a_fold01_intervals.png)

## Fold 2

![M2a fold 2 panel](figures/m2a_fold02_panel.png)

[Trace](figures/m2a_fold02_trace.png) · [Rank](figures/m2a_fold02_rank.png) · [Intervals](figures/m2a_fold02_intervals.png)

## Fold 3

![M2a fold 3 panel](figures/m2a_fold03_panel.png)

[Trace](figures/m2a_fold03_trace.png) · [Rank](figures/m2a_fold03_rank.png) · [Intervals](figures/m2a_fold03_intervals.png)

## Analysis role

Alternative fitted specification; not the model packaged in the M3b deployment bundle.

Each CV fold is a separate fit with its own chains. This batch does not yet include verified training/scoring entrypoints or their full dependency closure. No complete run command is claimed. PNG table export remains a separate reviewed addition.

## Individual fold pages

- [Fold 1: panel and individual plots](folds/fold01/README.md)
- [Fold 2: panel and individual plots](folds/fold02/README.md)
- [Fold 3: panel and individual plots](folds/fold03/README.md)
