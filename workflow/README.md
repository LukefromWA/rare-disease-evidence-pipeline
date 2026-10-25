# Analysis workflow

> This guide is shared between `main` and `workflow`. Run the commands from a checkout of the `workflow` branch; publishing this guide on `main` does not place the implementation there.

[Project overview](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/README.md) · [How to run](#how-to-run) · [CLI guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/CLI.md)

I use this workflow to carry case evidence through candidate ranking, frozen-model scoring and mechanism-aware drug-repurposing research. The calculations stay in the existing scripts and libraries; `bin/mva` provides the command-line entrypoint.

The workflow separates case analysis from external model training. The master workflow applies retained external-model artifacts without refitting those models on the case.

## Execution status

The configured Linux/WSL `full_extended` workflow completed a forced recomputation of all 62 scheduled jobs on 24 October 2026 in 27 minutes and 38 seconds, including the final completion and integrity audit.

The run covered core candidate analysis, frozen-model scoring, readthrough routing, transcript-context verification, genomic mapping, external-prediction lookups and source-linked claim ledgers.

The run used existing VEP annotations, prepared reference resources, retained external-model artifacts and authorized inputs. The recorded duration excludes installation, initial annotation, external-model training, resource acquisition and manual literature curation. It is an observed runtime for the configured environment, not a clean-install benchmark.

### Core and extended execution

The default `all` target requests resource readiness, A10 ranked candidates, the selected C4 handoff and the D7 report table. The `full_extended` target additionally connects the readthrough, external-prediction and claim-ledger stages to the final completion audit.

The CLI Track 1 target includes A11 in fresh-Stan mode. Standalone A11, optional Track 2 analyses and PDF reporting have separate execution paths. Completion of the core or extended target does not imply that PDF reporting was run.

### Earlier CLI validation: historical-handoff route

An earlier CLI validation completed Track 1, an explicitly prepared same-case historical candidate handoff and Track 2 using existing VEP output and prepared local resources. Output integrity and cross-stage consistency were checked, and all 15 CLI tests passed at that checkpoint.

| Stage | Completed steps | Observed elapsed time |
|---|---:|---:|
| Track 1 | 22/22 | 22 minutes, 7 seconds |
| Historical candidate handoff | 2/2 | About 5 seconds |
| Track 2 | 17/17 | 56 seconds |

These counts describe the earlier CLI route, not a breakdown of the later 62-job extended run.

The route was repaired and resumed during validation. These are observed single-core execution times, not a clean-install or uninterrupted runtime benchmark. Setup, resource preparation, raw annotation and manual handoff preparation are additional.

### Execution planning and reproduction

Before these analysis runs, the source-import checkpoint produced a successful 39-job Snakemake dry-run plan through A10, C4 and D7. That checkpoint verified the execution plan without running analysis jobs; the subsequent configured runs established execution completion.

Execution completion, environment restoration, integration-test coverage and clean-install reproduction are separate claims. Read the corresponding records alongside the configured run results.

## How to run

The instructions below cover the demonstrated CLI historical-handoff route. That route is distinct from the separately recorded master `full_extended` execution described above.

See the [execution status](#execution-status) for the completed core and extended run. The [CLI guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/CLI.md) documents path handling, configuration boundaries and the exact historical-handoff commands.

Use the route-specific configuration and documented entrypoint. Completion of one route does not establish completion of every optional analysis or reporting target.

### 1. Use the implementation branch and configured runtime

From an existing checkout with a clean working tree:

```bash
git fetch origin
git switch workflow
git merge --ff-only origin/workflow
```

Use Linux or WSL and activate your configured Python, R and Snakemake runtime. The launcher uses `python3` unless `MVA_PYTHON` names another interpreter. Snakemake and supporting tools must also be on `PATH`; selecting a Python interpreter alone does not install or activate the other dependencies.

```bash
# Optional: select the Python interpreter from your configured runtime.
export MVA_PYTHON=/path/to/runtime/bin/python
./bin/mva --help
```

Replace example paths with real local paths. This is not an automatic environment installer.

Provision the complete resource/model/analysis profile, frozen model artifacts, core resource bundle and matching Reactome hierarchy before execution. Keep private configurations, case inputs and generated results outside the source checkout.

### 2. Prepare and run Track 1

```bash
./bin/mva
```

Choose Track 1. Supply the complete profile, project/data root, cores, GRCh38 variant inputs, selected VCF sample, run label, observed HPO file and a new private workspace. Review the profile's QC thresholds, priors and mechanism assumptions before confirming them.

For the demonstrated existing-VEP route, provide the annotated `.vcf.gz`, matching `.tbi` index and real original VEP command/provenance file. The HPO file contains one `HP:nnnnnnn` term per line.

The wizard saves a private `run.yaml`, sample sheet and phenotype snapshots. It defaults to planning; execution requires explicit confirmation.

For a saved configuration, preview and then execute:

```bash
./bin/mva track1 --config /path/to/workspace/run.yaml --root /path/to/project-data --cores 1 --dry-run
./bin/mva track1 --config /path/to/workspace/run.yaml --root /path/to/project-data --cores 1
```

### 3. Prepare the historical handoff

Follow [Existing-case historical handoff preparation](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/CLI.md#existing-case-historical-handoff-preparation).

Verify the historical manifest, GRCh38 declarations, input checksums and same-case identity. Keep the completed Track 1 configuration unchanged, prepare a separate `track2_handoff.yaml`, and build a new union from the completed current-run artifacts using the documented restricted `layer_c` commands.

Track 1 alone does not prepare this historical union. The existing-run Track 2 wizard does not build it automatically. Do not substitute an older run's scored union for the new build. The current connected case adapter requires historical support.

### 4. Run Track 2

Use the prepared handoff configuration:

```bash
./bin/mva track2 --config /path/to/workspace/track2_handoff.yaml --root /path/to/project-data --cores 1 --dry-run
./bin/mva track2 --config /path/to/workspace/track2_handoff.yaml --root /path/to/project-data --cores 1
```

When continuing a saved case configuration, select `track1`, `track2` or `run` explicitly. Preserve the intended configuration and inspect the plan before continuing; do not treat continuation as a new independent validation run.

### 5. Independent candidate-gene exploration

This route does not require case VCF/HPO inputs or the historical handoff:

```bash
./bin/mva track2 --genes TP53 BRCA1 --resources /path/to/validated-core-resources --hierarchy /path/to/matching-reactome-hierarchy --out /path/to/new-gene-run --cores 1 --dry-run
```

A gene-only dry run saves its configuration, so its output directory then exists. Continue that saved run with:

```bash
./bin/mva resume --config /path/to/new-gene-run/run.yaml
```

Alternatively, omit `--dry-run` on the initial command and use a new output directory.

Gene-only discovery is unweighted direct-target/shared-pathway exploration, not case-based D0–D7 tiering. It does not create variant evidence, causal probabilities, Stan scores or mechanism weights.

### 6. Inspect outputs and run the CLI tests

Case outputs are under the configured output root and sample label. Start with:

- `a10/ranked_candidates.tsv.gz`: ranked candidates.
- `c3_stan_cv/stan_cv_predictions_summary.tsv`: frozen predictions.
- `a11_review/calculation/`: ranking-sensitivity results.
- `d3/*status.tsv`: method availability.
- `d5/held_or_excluded.tsv`: eligibility audit.
- `d6/status.tsv`: sensitivity and scenario coverage.
- `d7/final_report.tsv`: Track 2 report table.
- `logs/`: execution logs.

These paths describe the listed core and CLI outputs, not an exhaustive inventory of the extended stages. For the extended route, also inspect its completion/integrity audit and the recorded readthrough, transcript-context, mapping, external-prediction and claim-ledger outputs.

Gene-only outputs include `report/gene_summary.tsv` and `report/report.json` in that route's output directory.

A case report with no eligible drugs can be a valid result. Inspect the hold/exclusion reasons and method coverage rather than assuming every method ran or every candidate passed the gates.

From the repository root, use the configured Python:

```bash
"${MVA_PYTHON:-python3}" -m unittest discover -s tests -p 'test_cli*.py' -v
```

## Entrypoints and targets

| Entrypoint | Role |
|---|---|
| `bin/mva` | Interactive and noninteractive CLI, run from the repository root |
| `workflow/Snakefile` | Configured no-training master workflow |
| `workflow/track2_genes.Snakefile` | Independent gene-only discovery |
| `workflow/layer_a/Snakefile` | Standalone A11 sensitivity/reporting |

Rules live in `workflow/rules/`; calculations live in `workflow/scripts/` and `workflow/lib/`. Layer directories are not six interchangeable independent runners.

### Master targets

| Target | Scope |
|---|---|
| `track1` | A10, frozen Stan predictions and A11 in fresh-Stan mode |
| `layer_a` | Through A10 |
| `layer_c` | Selected C4 handoff |
| `layer_d` | Through D7 |
| `resources` | Native preparation and core resource validation |
| `all` | Resource readiness, A10, selected C4 handoff and D7 |
| `full_extended` | Core route plus connected readthrough, external-prediction and claim-ledger stages, ending with the completion and integrity audit |

`layer_b` requests an existing artifact; it does not train a model. The master does not import external-training rules. Fresh-Stan scoring mode does not mean that the external models are refitted on the case.

E21 writes a separate validation-status report and does not imply method enablement. Resource readiness, method availability and completed analysis are distinct statuses.

## Analysis boundaries

### Case evidence and frozen inference

Case evidence preparation, external model development and frozen inference remain separate. The case is not used to refit the external Stan models.

Candidate registries preserve their evidence sources and availability. Genomic evidence, phenotype similarity, network scores, sensitivity frequencies and model probabilities are not interchangeable quantities.

### Track 2 research

Track 2 examines mechanism scenarios, targets, compound evidence, eligibility gates and ranking sensitivity. A completed analysis does not establish biological causality, drug efficacy, pediatric safety or clinical suitability.

Method-status records distinguish available, disabled, unavailable and not-applicable analyses. Missing evidence is not a successful check or a negative biological result.

### Readthrough and external predictions

The connected readthrough extension adds transcript-context verification, genomic mapping, published RTDetective prediction lookups and source-linked literature evidence.

Computational predictions, transcript-context checks and literature claims remain distinguishable from demonstrated biological effects. The claim ledger records source-linked evidence rather than converting a prediction into an experimental result.

## Validation and execution records

The configured Linux/WSL workflow completed all 62 scheduled jobs in the
`full_extended` run, including the final completion and integrity audit.
The run used existing VEP annotations, prepared reference resources,
retained model artifacts and authorized inputs.

Additional validation records include:

- 15 passing CLI tests at the earlier validation checkpoint.
- Execution of the independent gene-only discovery route.
- A regression check confirming that all five weighted D2 discovery
  tables remained content-identical.

These checks cover different execution paths and are recorded separately.
Method-status outputs, evidence coverage and eligibility audits document
which analyses contributed to each result.

### Execution scope

The completed run begins with existing VEP annotations. Raw VCF-to-VEP
execution and the full interactive route remain separate integration-test
targets.

D7 produces the research report table. Optional PDF reporting follows a
separate execution path; its legacy asset wiring remains to be aligned
with the fresh-Stan handoff.

## Data and reproduction requirements

Run the workflow with the documented runtime, matching reference
resources, retained model artifacts, configuration and authorized inputs.
The execution records describe the tested environment; fresh-checkout
setup is documented separately.

Private case inputs, populated private configurations and generated case
results are not distributed with the public source.

## Intended use

This workflow supports computational research and human review.
Candidate rankings, model predictions and drug-repurposing hypotheses
