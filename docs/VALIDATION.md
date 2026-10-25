# Validation and reproducibility

Validation combines synthetic regression checks, command-line tests, and
configured workflow execution. Each addresses a different part of the system.

## Configured workflow execution

The configured Linux/WSL `full_extended` workflow completed a forced
recomputation of 62 scheduled jobs in 27 minutes and 38 seconds, including
the final completion and integrity audit.

The run used existing VEP annotations, prepared reference inputs, retained
external-model artifacts, local prediction tracks, and saved curated-source
bundles. The recorded duration excludes initial annotation, external-model
training, resource acquisition, and manual literature curation.

The completion audit checks declared output integrity, identity consistency,
coverage, provenance, and completion requirements.

## CLI and integration checks

Earlier validation exercised existing-case CLI dispatch, independent gene-only
discovery, and regression checks showing content-identical weighted D2 outputs.
That configured route was repaired and resumed during validation.

Public CLI regression tests include:

- `tests/test_cli_resource_manifest.py`
- `tests/test_cli_interactive.py`

The configured backend execution and command-line tests have distinct scopes.
Raw VCF-to-VEP execution, every interactive route, and portable fresh-checkout
reproduction require their own integration evidence.

## A11 regression coverage

Synthetic checks covered:

- Original calculation and non-label report-value regression.
- Fixed-seed repeatability and candidate coverage preservation.
- Input/provenance mismatch and overwrite rejection.
- Separate Stan and A11 source weights.
- Coverage categories, neutral rank labels, and CLI edge cases.
- Standalone calculation, reporting, and union attachment.
- A no-work second dry run.

These checks used reduced simulation counts without changing production
defaults. The internal A11 regression scripts and execution logs are retained
locally; they are separate from the public CLI tests.

## Interpretation

Software completion and numerical regression do not establish biological
accuracy, clinical suitability, or treatment efficacy. Synthetic regression
does not demonstrate real-case predictive performance.

Recorded execution demonstrates the configured route, not automatic resource
acquisition or a clean-install benchmark. Reproduction requires matching
inputs, resource releases, model artifacts, configuration, and runtime.

## Further reading

- [CLI usage and execution record](CLI.md)
- [Extended workflow and audit](full_extended.md)
- [A11 execution and interpretation](A11_REPRODUCIBILITY.md)
