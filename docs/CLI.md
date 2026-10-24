# CLI usage and tested scope

Run `bin/mva` without arguments for the interactive menu.
Use the tested Linux/WSL Python, R and Snakemake environment on PATH.

## Interactive choices

1. Track 1: case evidence, frozen scoring and A11.
2. Track 2: existing run or manually supplied candidate genes.
3. Connected Track 1 to Track 2 workflow.

New-case setup asks whether raw VCF and existing VEP output are available.
Existing output requires an indexed annotated `.vcf.gz` and its real VEP
command/provenance file. Raw mode requires an indexed VCF, GRCh38 FASTA/index,
and an installed VEP cache/version. No provenance is fabricated.

Track 1 requires a complete resource/model/analysis profile. Its QC thresholds,
priors and resources are retained after explicit user confirmation. This is not
an automatic environment installer or a substitute for reviewing those settings.

The wizard checks VCF sample names, asks for the pipeline sample label and observed
HPO terms, replaces old private input fields, and writes a private `run.yaml`,
sample sheet and phenotype snapshots outside the source repository.
It defaults to planning rather than execution.

## Paths

Linux paths and Windows drive paths are accepted under WSL, including paths
copied with surrounding quotes. Windows paths are translated using `wslpath`
when available. Native Windows workflow execution is not validated.
Mount network shares first and use their Linux paths.
Control characters and commas in paths are rejected because of legacy argument
contracts. On the command line, quote Windows paths so Bash preserves backslashes.
Interactive pasted paths do not pass through Bash parsing.

## Noninteractive commands

- `bin/mva track1 --config /path/run.yaml --root /path/data --dry-run`
- `bin/mva run --config /path/run.yaml --root /path/data --dry-run`
- `bin/mva track2 --from-run /path/case-workspace --root /path/data --dry-run`
- `bin/mva track2 --genes TP53 BRCA1 --resources /path/core --hierarchy /path/hierarchy --out /path/new-run --dry-run`
- `bin/mva resume --config /path/gene-run/run.yaml`

A gene-only dry run creates its saved configuration but does not execute discovery.
For master runs, select track1/track2/run explicitly when continuing a saved config.
The generated master config records the project root, so it need not be retyped.

## Boundaries

Track 1 does not require the historical union: it requests A10, frozen Stan
predictions and A11 directly. Connected mode still needs explicitly declared,
same-case historical support for the current union adapter. A no-history connected
adapter is not implemented; the wizard offers Track 1 alone instead.

Gene-only mode reuses D2 direct-target and shared-pathway discovery. It creates
no variant coordinates, causal probabilities, Stan scores or mechanism weights.
Unmapped inputs remain explicit. Its report is unweighted exploration, not D0-D7
research tiering. Mechanism-specific imports, label review and safety adjudication
are not implemented in this initial gene-only route.

Family-history terms are retained separately; this does not imply they enter
the current proband-only phenotype calculation.

## Validation record

The configured backend completed after repair/resume, with integrity and
cross-stage checks. Existing-case CLI dispatch, gene-only discovery, and
content-identical weighted D2 regression were demonstrated locally.
Manifest-reader regression tests are in tests/test_cli_resource_manifest.py.
Interactive configuration/path tests are in tests/test_cli_interactive.py.
Raw VCF-to-VEP execution and full interactive end-to-end execution still require
dedicated integration tests. A syntax/configuration test is not that validation.

## Repository organization

`workflow/Snakefile` is the canonical configured no-training master.
`workflow/track2_genes.Snakefile` is the independent unweighted gene entrypoint.
`workflow/layer_a/Snakefile` is standalone A11, not another complete Layer A runner.
The scientific script/library locations remain stable.

Generated configs, gene lists, model objects, payloads and case reports must not
be committed. No training, automatic safety approval or clinical recommendation
is introduced by the CLI.

## Existing-case historical handoff preparation

The tested route runs Track 1, prepares the historical candidate union, then
runs Track 2. Observed execution times were about 22 minutes, 5 seconds and
1 minute respectively, using one core and prepared local resources.

### 1. Verify the historical source

Use historical C4/A10 files from the same case. The historical union directory
must contain `manifest.json` with:

- Version `A10_C4_union_Stan_v1_2`.
- `current_build` and `legacy_build` both set to `GRCh38`.
- Exactly one `legacy_c4` and one `legacy_a10` input.
- SHA-256 values matching those input files.

Confirm case identity separately from checksum verification.

### 2. Save a handoff configuration

Keep the completed Track 1 `run.yaml` unchanged. Copy its full configuration
to a private `track2_handoff.yaml` and set the sample-specific historical source:

```yaml
track2:
  frozen_union_dirs:
    SAMPLE_LABEL: /path/to/verified-historical-union
```

Merge this setting into the existing `track2` section; the fragment above is
not a complete configuration. Keep the sample label, output root and other
analysis settings unchanged. Verify that only `frozen_union_dirs` differs.

### 3. Build the new union from completed Track 1 outputs

Activate the configured runtime, then set:

```bash
export MVA_ROOT=/path/to/project-data
export MVA_CONFIG=/path/to/workspace/track2_handoff.yaml
export MVA_REPO=/path/to/repository
```

Preview the handoff:

```bash
snakemake \
  --snakefile "$MVA_REPO/workflow/Snakefile" \
  --directory "$MVA_ROOT" \
  --cores 1 --printshellcmds --dry-run \
  layer_c --allowed-rules c4_fresh_stan_union layer_c
```

The plan should contain only `c4_fresh_stan_union` and `layer_c`.
The restricted build explicitly reuses completed Track 1 artifacts.
Check that their configuration and provenance match the intended analysis.

Execute the reviewed plan:

```bash
snakemake \
  --snakefile "$MVA_REPO/workflow/Snakefile" \
  --directory "$MVA_ROOT" \
  --cores 1 --printshellcmds --rerun-incomplete \
  layer_c --allowed-rules c4_fresh_stan_union layer_c
```

This builds a new union using the current run's rankings and predictions.
Do not copy an older run's scored union into the new results directory.

### 4. Run Track 2

```bash
"$MVA_REPO/bin/mva" track2 \
  --config "$MVA_CONFIG" \
  --root "$MVA_ROOT" \
  --cores 1
```

This route checks the handoff and resources before executing Layer D.
The existing-run wizard does not automatically perform the preparation above.

### 5. Check the report and coverage

Inspect `d7/final_report.tsv`, `d3/*status.tsv`,
`d5/held_or_excluded.tsv` and `d6/status.tsv` under the sample's output
directory. Optional methods record their availability. A report with no
eligible drugs retains the reasons for holds and exclusions.

The current connected adapter requires historical support. Independent
candidate-gene exploration is available without that handoff, but is unweighted
and is not equivalent to the case-based D0-D7 analysis.

Keep private configurations, historical case artifacts and generated results
outside the public repository.
