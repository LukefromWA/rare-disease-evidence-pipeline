# A/B/C source checkpoint — 2026-10-22

## Scope
Selected working sources copied unchanged. Existing A core, Stan definitions,
documentation and tests remain unchanged. This is a source checkpoint,
not a complete A0–A11 or B/C end-to-end execution release.
Layer D is deferred and is not a prerequisite for this checkpoint.

## Statistical paths
- b11_fit_m2.R, b14_freeze.R, c1_features.R and c3_frozen_inference.R
  include the separate optimization-model baseline interface.
- b14_stan_cv_pipeline.R and b14_stan_cv_helpers.R package and score
  an equal-weight mixture of three retained M3b fold components.
  This is not a full-data Bayesian posterior.
- c1_evidence_handoff_v1.r preserves evidence without changing predictions.
- c4_d0_union_stan_v1_2.R retains current A10 UNION identifiable historical
  C4. Historical-only candidates remain protected and unweighted.
  Unresolved historical rows remain separately documented.

## Remaining integration work
- Verify current G/H/N definitions against original training inputs.
- Verify genuine four-chain metadata in retained artifacts.
- Correct or explicitly override default source/helper/Stan paths.
- Preserve divergence and unverified-treedepth qualifications.
- Harden output preflight and rerun behavior.
- Correct older runner diagnostics and score_frozen API usage.
- Record effective scientific configuration without private overlays.
- Add layer/master workflow entrypoints and connected tests.
- Keep normal case analysis separate from explicitly requested training.
- Carry protected-unscored candidate coverage into downstream processing.
- Existing baseline c1_features.R has a fixed case label; it is not
  the general case-input wrapper.

## Boundaries
No populated configurations, patient inputs, phenotype lists, fitted chains,
predictor bundles, predictions, case tables, logs or runtime manifests
are imported. The complete A runner and A11 are not newly imported here.
Original working sources are not modified.
No fitting, scoring, calibration or new tests were run by this import.
Staging does not constitute new privacy certification.
Existing test evidence remains in LOCAL_TEST_RESULTS.md.

## Input contract
Proband HPO text: one HP identifier per line. Family history separate.

## Git policy
The repository deny-by-default ignore policy is unchanged.
Only explicitly listed source and checkpoint files are force-added.
