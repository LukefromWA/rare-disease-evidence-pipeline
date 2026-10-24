# Candidate-scoped streaming A1

A0 -> a0_candidate_subset -> a1_streaming_selftest -> a1_evidence_matrix -> A2.
A3 remains unchanged in this change set. No VEP rerun is introduced.

## Policy

The default `legacy_any_csq_token_position_v1` deliberately reproduces the legacy
algorithm: match an entire HIGH/MODERATE token in any CSQ field, case-insensitively,
then retain every record at selected CHROM/POS sites. This is not schema-aware
IMPACT-only filtering, nor exact-allele filtering. A non-IMPACT field containing
HIGH and other alleles at a selected site can affect membership. Those behaviors
are covered by synthetic regression tests, not silently corrected.

All surviving INFO, FORMAT and CSQ annotations are retained. No AF, DP, GQ, AD,
VAF, FILTER, QUAL or contig exclusion is added upstream. Candidate rows failing
core QC remain in A1. Counting is accumulated in the selection pass. Progress
is logged every configured number of scanned records. Extraction still makes
another full-input pass. No performance multiplier is guaranteed.

## Optional configuration

Merge `config/a1.candidates.example.yaml` into a complete run configuration; the
fragment is not a standalone workflow configuration. Candidate policy, batch size
and progress interval have explicit defaults in the rule module. Existing QC,
AF, representative-transcript and sample-selection configuration is unchanged.

## Outputs and migration

Candidate input: `<output_root>/<sample>/a0_candidates/` contains candidates.vcf.gz,
its tabix index, and candidate_manifest.json. The manifest reports actual policy,
counts, hashes, tool version and the absence of upstream hard-QC/contig filtering.
A1 retains its standard a1/ output path and adds candidate_scope.json, a copy of
the producer manifest. The A1 manifest independently fingerprints its VCF input.
The copied scope file is a companion provenance artifact, not a historical-parity
certificate. Existing unmarked A1 products block DAG construction when A1 is
requested: preserve them and choose a fresh run.output_root. Do not force-delete
old full-WGS matrices to migrate. Downstream products also belong in the fresh
output root. Incomplete first-time candidate A1 products require explicit review.

The builder is bounded by record count, not a strict byte-level RAM guarantee.
External ordering uses a 128 MiB sort budget. Long CSQ records can still be large.
The current builder stops on empty candidate VCFs; the selector gives an explicit
no-candidates diagnostic and publishes no final products in that case.

## Validation

Run `python3 tests/test_candidate_subset.py` with bcftools, awk and gzip installed.
These synthetic tests cover case folding, token boundaries, annotations on
multiple transcripts, position-level collateral retention, INFO/FORMAT/phase
retention, low-quality-row retention, refusal to overwrite, and zero candidates.

The DAG automatically self-tests the active R builder on a bounded prefix at
two batch sizes and checks duplicate rejection before building A1. It is not
proof of historical full-dataset equivalence. R helper/code/config hashes bind
the parity report; it must be regenerated after code or configuration changes.
Self-test runs use a new child directory and publish only the passed report to
the declared output, avoiding Snakemake's automatic output-parent creation.

A private validation already demonstrated complete allele-scope recovery from
an annotated full VCF and a successful compact build. No patient inputs, private
paths, sample identifiers, case counts or saved reference VCFs are repository
fixtures or required production dependencies. New source locations generate
new manifests and parity hashes. Use the pinned R environment for reproduction.

The current A3 empty-HGVS overwrite behavior is intentionally deferred to a
separate reviewed change. No A3 semantic changes are bundled here.
