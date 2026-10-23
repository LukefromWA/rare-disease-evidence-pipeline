# A11 reproducible execution path

This port preserves the original 18-specification Dirichlet/Beta sensitivity
calculation: three scenarios, two network views and three gamma values.
Production defaults remain 10,000 draws per specification and seed 20260912.
This is not a calibrated causal/pathogenicity probability model.
Dependence is used for network-weight adjustment; the original Beta score
perturbations are not newly replaced with joint copula sampling.

## Ownership
- A: workflow/layer_a/scripts/a11_dependence_sensitivity.R
- C: workflow/layer_c/scripts/attach_a11_multimodel.R
- F: workflow/layer_f/scripts/f_a11_coverage_report.R
- Shared execution/provenance: workflow/lib/a11_runtime.R
- Standalone runner: workflow/layer_a/Snakefile

The standalone runner starts from an existing candidate-evidence CSV.
It does not implement the missing upstream A10-to-evidence builder or the
full A0-F master. An existing validated union is needed for attachment.

## Input
Calculation requires candidate_id, gene_symbol, genomic, hpo, string,
multiplex. Report additionally requires has_genomic, has_hpo, has_string,
has_multiplex. Candidate IDs must be unique/nonblank, gene symbols nonblank;
at least two candidates must survive the unchanged two-view eligibility
filter. That minimum fails closed on an unsupported original shape edge case.
Input row order remains unchanged because it affects Monte Carlo draws.

## Execution
Use config/a11.example.yaml as a template, keeping populated paths private:

    snakemake -s workflow/layer_a/Snakefile --configfile PRIVATE_CONFIG --cores 1

Without attachment enabled, only calculation and coverage reporting run.
Direct Rscript entrypoints accept explicit --input/--outdir paths; the
calculation also accepts --seed, --draws and optional --focus-gene.
The report requires --consensus, --copula and --a11-manifest.
Attachment retains --union, --consensus, --a11_input and --outdir and now
requires --a11_manifest. No runtime package installation is performed.

## Preserved and changed behavior
Formulas, priors, scenario order, candidate filtering, correlation fallback,
shrinkage, missing-score handling, quantiles and ranking are preserved.
Legacy posterior_* column names are retained for compatibility, not causal
interpretation. Consensus is an equal arithmetic average over sensitivity
specifications, not learned model probabilities.
Fixed gene-focused exports are generalized to optional --focus-gene reporting;
no gene is selected by default. Full candidate outputs remain unchanged.
The existing top-5/top-10/top-25 exports are review conveniences, not the
complete candidate registry or a downstream exclusion rule.

Each producer/report run writes run_provenance.json with source/helper/input
hashes, settings, package versions, R version, RNG kind, locale and output
hashes. Producer runs also save initial/final RNG states and an eligibility
ledger. The new report/attachment verify hashes against producer provenance.
Existing legacy artifacts require honest provenance reconstruction or a
separate legacy path; do not fabricate provenance for them.

Fresh output paths are required by the R entrypoints. A failed run may leave
an incomplete directory; absence of final provenance means incomplete.
The Snakemake outputs are dedicated directory outputs, never shared with
other rules. Runtime outputs can be patient-derived and must remain outside
public source exports. No populated config or real case result is committed.

## Validation and limitations
See A11_REWIRE_TEST_RESULTS.md for tests actually run by the patcher.
Synthetic regression tests compare original and ported calculations/reports
within one installed environment; they do not validate biology or prove
cross-platform bitwise equality. Timestamps are not repeatability criteria.
A repository-wide locked environment and the full master DAG remain pending.

## Original source hashes
- Bayesian_Dependence_Aware_Variant_Ranking_final.r: `d7d4c54db0073a59af0d3d70f585054e373ad0bb7b1a48d96e72f21cb4b166d8`
- 13_make_bayesian_coverage_report.r: `d06cf6cd501c35ecd180f4e27ac6ff5c3c719388358dedc1a08872214b6ebf38`
- attach_a11_multimodel_v1_1.R: `33070452d8bee77c741650422b4f547cc0e540e49cdd8996fb7a475938185102`
