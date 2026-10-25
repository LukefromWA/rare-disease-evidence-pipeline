# Analysis workflow

> This guide is shared between `main` and `workflow`. Run the commands from a checkout of the `workflow` branch; publishing this guide on `main` does not place the implementation there.

[Project overview](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/README.md) · [How to run](#how-to-run) · [CLI guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/CLI.md)

I use this workflow to carry case evidence through candidate ranking, frozen model scoring and mechanism-aware drug-repurposing research. The calculations stay in the existing scripts and libraries; `bin/mva` provides the command-line entrypoint.

## Earlier CLI validation: historical-handoff route

The configured Linux/WSL route completed Track 1, an explicitly prepared same-case historical candidate handoff, and Track 2 using existing VEP output and prepared local resources. Output integrity and cross-stage consistency were checked, and all 15 CLI tests passed.

| Stage | Completed steps | Observed elapsed time |
|---|---:|---:|
| Track 1 | 22/22 | 22 minutes, 7 seconds |
| Historical handoff | 2/2 | About 5 seconds |
| Track 2 | 17/17 | 56 seconds |

The workflow was repaired and resumed during validation. These are observed single-core execution times, not a clean-install or uninterrupted runtime benchmark. Setup, resource preparation, raw annotation and manual handoff preparation are additional.

## How to run

The instructions below cover the demonstrated historical-handoff route. See the [project execution status](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/README.md#current-execution-status) for the separately recorded core and extended workflow execution. The [CLI guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/CLI.md) includes path handling, configuration boundaries and the exact historical-handoff commands.

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

Replace example paths with real local paths. This is not an automatic environment installer. Provision the complete resource/model/analysis profile, frozen model artifacts, core resource bundle and matching Reactome hierarchy before execution. Keep private configurations, case inputs and generated results outside the source checkout.

### 2. Prepare and run Track 1

```bash
./bin/mva
```

Choose Track 1. Supply the complete profile, project/data root, cores, GRCh38 variant inputs, selected VCF sample, run label, observed HPO file and a new private workspace. Review the profile's QC thresholds, priors and mechanism assumptions before confirming them.

For the demonstrated existing-VEP route, provide the annotated `.vcf.gz`, matching `.tbi` index and real original VEP command/provenance file. The HPO file contains one `HP:nnnnnnn` term per line. The wizard saves a private `run.yaml`, sample sheet and phenotype snapshots. It defaults to planning; execution requires explicit confirmation.

For a saved configuration, preview and then execute:

```bash
./bin/mva track1 --config /path/to/workspace/run.yaml --root /path/to/project-data --cores 1 --dry-run
./bin/mva track1 --config /path/to/workspace/run.yaml --root /path/to/project-data --cores 1
```

### 3. Prepare the historical handoff

Follow [Existing-case historical handoff preparation](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/CLI.md#existing-case-historical-handoff-preparation). Verify the historical manifest, GRCh38 declarations, input checksums and same-case identity. Keep the completed Track 1 configuration unchanged, prepare a separate `track2_handoff.yaml`, and build a new union from the completed current-run artifacts using the documented restricted `layer_c` commands.

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

Alternatively, omit `--dry-run` on the initial command and use a new output directory. Gene-only discovery is unweighted direct-target/shared-pathway exploration, not case-based D0-D7 tiering. It does not create variant evidence, causal probabilities, Stan scores or mechanism weights.

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

Gene-only outputs include `report/gene_summary.tsv` and `report/report.json` in that route's output directory. A case report with no eligible drugs can be a valid result; inspect the hold/exclusion reasons and method coverage rather than assuming every method ran.

From the repository root, use the configured Python:

```bash
"${MVA_PYTHON:-python3}" -m unittest discover -s tests -p 'test_cli*.py' -v
```

## Entrypoints and targets

- `bin/mva`, from the repository root: interactive and noninteractive CLI.
- `workflow/Snakefile`: configured no-training master.
- `workflow/track2_genes.Snakefile`: independent gene-only discovery.
- `workflow/layer_a/Snakefile`: standalone A11 sensitivity/reporting only.

Rules live in `workflow/rules/`; calculations live in `workflow/scripts/` and `workflow/lib/`. Layer directories are not six interchangeable independent runners.

Master targets:

- `track1`: A10, frozen Stan predictions and A11 in fresh-Stan mode.
- `layer_a`: through A10.
- `layer_c`: selected C4 handoff.
- `layer_d`: through D7.
- `resources`: native preparation and core resource validation.
- `all`: resource readiness, A10, selected C4 handoff and D7.

`layer_b` requests an existing artifact; it does not train a model. The master does not import external-training rules. E21 writes a separate validation-status report and does not imply method enablement.

## Validation boundary

The configured execution does not establish portable fresh-checkout reproduction or clinical validity. Raw VCF-to-VEP execution and full interactive end-to-end execution still need dedicated integration tests. The CLI gene-only route was executed, and all five weighted D2 discovery tables remained content-identical in a regression run.

D7's safety-gated research table is separate from optional PDF reporting. Legacy PDF asset wiring still needs alignment with the fresh-Stan handoff; PDF generation is not part of the demonstrated `all` target. Outputs support research and human review, not treatment recommendations.
