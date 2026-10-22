# Rare-Disease Evidence Pipeline

Private development repository for evidence-preserving rare-disease
candidate prioritization and mechanism-aware research exploration.

## Status

Source consolidation and reproducibility review are in progress.
This initial repository contains documentation and ignore rules only.
Reviewed code and synthetic tests will be imported in separate commits.

## Logical layers

- A: private-case evidence processing and candidate ranking.
- A11: ranking sensitivity under declared assumptions.
- B: external-cohort preparation and model development.
- C: frozen model application without private-case refitting.
- C4/D0: candidate registry and mechanism-scenario handoff.
- D: mechanism-aware discovery, gating and sensitivity.
- E: case-independent resource preparation and validation.
- F: run-specific reporting and provenance.
- Release: reviewed export and packaging controls.

## Interpretation

Normalized evidence weights, A11 rank frequencies and Stan predictive
probabilities are different quantities and must remain separately labelled.

The inspected Stan deployment adapter uses an equal-weight mixture of
three retained cross-validation predictors, not a single full-data posterior.
Provisional sampling qualifications must remain visible.

Unavailable or disabled analyses must not be described as completed.

## Data boundaries

Do not upload patient inputs, phenotype lists, private configuration,
case-derived results, prediction objects, candidate registries, generated
reports, credentials or the private consolidation workspace.

Repository privacy does not authorize uploading protected data.
Shared development should use reviewed code, schemas and synthetic fixtures.

## Validation

Three export-boundary tests passed in the local development workspace.
This is not complete privacy clearance or end-to-end scientific validation.
The tests and reviewed implementation have not yet been imported here.

## Intended use

Computational research and human review only.
Not a clinically validated diagnostic or treatment-selection system.

## License

Not yet selected. Third-party materials retain their own licenses.

## AI assistance

Perplexity (Pro tier) was used primarily for documentation, debugging and
code/workflow integration. All relevant privacy
and data-usage settings preventing private data from being used for model
training were enabled from the beginning and checked throughout the project.
Equivalent training-use restrictions.

See [AI assistance disclosure](docs/AI_ASSISTANCE.md). Provider/product/tier
details for any additional contributing service belong in the final methods
record. This statement does not assert zero retention.

