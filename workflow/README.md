# Workflow

[Project overview](../README.md) · [Model specifications](../models/README.md) · [Workflow implementation](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow)

The integrated workflow is maintained on the `workflow` branch. This directory
on `main` provides navigation and a summary of the analysis layers; execution
instructions and source code are maintained alongside the implementation.

## Implementation and execution guides

| Looking for | Start with |
|---|---|
| Setup, analysis stages, and step-by-step execution | [Workflow guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/workflow/README.md) |
| CLI inputs, configuration, and execution routes | [CLI guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/CLI.md) |
| Master workflow and targets | [Snakefile](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/workflow/Snakefile) |
| Stage dependencies and commands | [Rule modules](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow/rules) |
| Executable analysis stages | [Scripts](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow/scripts) |
| Shared calculations and validation | [Workflow library](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow/lib) |
| Readthrough, external predictions, and source-linked claim ledgers | [Extended workflow guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/full_extended.md) |
| External-model specifications and retained-fit records | [External-training branch](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/external-training) |

## Execution status

The configured Linux/WSL `full_extended` workflow completed a forced
recomputation of all 62 scheduled jobs on 24 October 2026 in 27 minutes
and 38 seconds, including the final completion and integrity audit.

The run covered core candidate analysis, frozen-model scoring, readthrough
routing, transcript-context verification, genomic mapping, external-prediction
lookups, and source-linked claim ledgers.

Execution uses the documented runtime, prepared reference resources, retained
model artifacts, and authorized inputs. See the project overview and execution
guides for detailed run conditions, earlier CLI validation, and method-specific
availability.

## Layer responsibilities

| Layer | Role and handoff |
|---|---|
| A | Prepare case evidence, apply QC, compare phenotype/network views, integrate evidence, and retain candidate handoffs |
| A11 | Examine ranking sensitivity under declared perturbations, separately from causal probability |
| B | Prepare external public cohorts, audit deduplication and splits, fit features, develop models, and freeze retained artifacts |
| C | Apply the frozen feature/model contract without case-specific refitting |
| A ∪ C | Reconcile candidate identities and preserve both source views and their distinct score meanings for Track 2 |
| D | Examine mechanism scenarios, complementary discovery evidence, evidence integration, gates, and sensitivity |
| E | Prepare, version, and validate case-independent reference resources |
| F | Assemble run-specific tables, figures, provenance, and reports |

## Analysis design

External-model development remains separate from case execution. The case
workflow applies retained model artifacts without refitting them on the
challenge genome.

The external-training sequence preserves B8 → B10 → B9: candidate construction,
group-aware split assignment, and training-partition feature fitting. Frozen
inference follows the matching feature and preprocessing contract.

Track 2 receives the union of retained A and C candidates while preserving
candidate identities, source membership, evidence coverage, and separate score
meanings. Ranking sensitivity, external predictions, and mechanism-led research
remain distinct from demonstrated biological effects or treatment efficacy.
