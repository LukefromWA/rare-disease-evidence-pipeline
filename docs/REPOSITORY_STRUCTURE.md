# Repository organization

The repository separates the project overview, integrated implementation, and
external-model development.

## Branches

| Branch | Purpose |
|---|---|
| `main` | Project overview, model descriptions, acknowledgements, and essential project records |
| `workflow` | Integrated analysis source, command-line tools, execution guides, and reference configuration |
| `external-training` | External-model source, retained-fit documentation, and diagnostic material |

The A–F layers describe analysis responsibilities within the integrated workflow;
they are not separate development branches.

## Implementation layout

| Location on `workflow` | Responsibility |
|---|---|
| `bin/mva` | Interactive and noninteractive command-line entrypoint |
| `config/` | Reference configuration and parameter records |
| `environment/` | Dependency and environment records |
| `tests/` | Public regression tests |
| `workflow/Snakefile` | Configured no-training master workflow |
| `workflow/track2_genes.Snakefile` | Independent gene-only exploration |
| `workflow/rules/` | Stage dependencies and execution commands |
| `workflow/scripts/` | Executable analysis stages |
| `workflow/lib/` | Shared calculations, validation, and provenance helpers |
| `workflow/stan/` | Bayesian model definitions |
| `workflow/layer_a/` through `workflow/layer_f/` | Layer-specific components and entrypoints where provided |
| `models/` | Model specifications and interpretation |
| `docs/` | Execution guides, methods, citations, and disclosures |

The top-level `lib/` directory is separate from the implemented workflow library
in `workflow/lib/`. Layer directories do not each provide a complete independent
pipeline.

## Analysis responsibilities

- A prepares case evidence, ranks candidates, and examines sensitivity.
- B develops models from external public cohorts.
- C applies frozen model artifacts and reconciles candidate views.
- D explores mechanisms, targets, compound evidence, and research eligibility.
- E prepares and validates reference resources.
- F assembles saved results and provenance for review.

External training is separate from case execution. The case workflow applies
retained model artifacts without refitting them on the case.

## Data boundaries

Authorized case inputs, populated private configurations, generated case results,
and local resource/model bundles remain outside the public source checkout.

## Navigation

- [CLI guide](CLI.md)
- [Extended workflow](full_extended.md)
- [Standalone A11 analysis](A11_REPRODUCIBILITY.md)
- [Validation scope](VALIDATION.md)
