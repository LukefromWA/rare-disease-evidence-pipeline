# Full extended workflow

The explicit `full_extended` target includes the configured core outputs, readthrough routing,
CDS context verification, genomic mapping, RTDetective coordinate lookup/adoption,
source-anchored claim ledger, mechanistic augmentation and per-sample completion audit.
The default `all` target remains the core workflow.

## Configuration and execution
Keep the normal case configuration as MVA_CONFIG and provide a typed YAML overlay based on
config/full_extended.example.yaml. Replace its sample/gene/path placeholders. A YAML overlay
avoids nested command-line Boolean string conversion. Supply your resource files separately;
no case outputs, source snapshots, large tracks or machine-specific overlays are distributed.

```bash
export MVA_ROOT=/path/to/project-data
export MVA_CONFIG=/path/to/case/run.yaml
snakemake --snakefile workflow/Snakefile --directory "$MVA_ROOT" --cores 1 \
  --printshellcmds --rerun-incomplete full_extended \
  --configfile /path/to/extended.overlay.yaml --dry-run
```
Review the plan, then repeat without --dry-run. Do not suppress provenance rerun checks.
The requested sample must exist in the base samples file and the selected genes must be
compatible with the context producer. Context scope should reflect reviewed routing;
it is not an automatic general-purpose selector. The aggregate currently expects fresh
c4_stan_union outputs and the reviewed three-track adoption profile.

## Input contracts
- Fresh C4 variant evidence and A10 ranked candidates feed the codon context producer.
- GRCh38 positive-strand coding SNVs are supported by the mapping/lookup adapter. Flank
  windows may cross exons, but the stop codon must map to three contiguous genomic bases.
- Mapping validates the full genomic reference window against the frozen CDS and the exact
  candidate alleles. Transcript-version checks do not independently guarantee immutable
  version binding by the mapping endpoint; the manifest retains this limitation.
- The three RTDetective UGA tracks require the pinned inventory for article 23708901 v2.
  The scanner verifies published size/MD5, records SHA256, preserves raw score strings and
  source-line numbers, and rejects malformed/nonfinite rows. It does not reverse descending
  intervals or aggregate duplicates. Nonforward overlap blocks unique-match adoption.
- Saved fulltext and bibliography bundles require output_sha256 manifests. Claim seeds use
  schema_version 1 and source anchors checked by the ledger. The mechanistic bundle uses
  primary_mechanistic_source_collection_v1. These are explicit curated inputs, not a claim
  of independent compound discovery.
- Reference retrieval is opt-in. It can disclose transcript IDs and genomic coordinates to
  external services. mechanism-search offline mode does not control these reference flags.
  The track scanner, ledgers and completion audit themselves are local-only.

## Outputs and verification
The completion artifact is <output_root>/<sample>/extended/run_audit.json. It checks bundle
hashes, routing/bridge identity coverage, three-track adoption statuses, prediction/status
agreement, lookup mapping fingerprints and claim-review coverage. No-hit/review statuses are
reported rather than replaced with predictions. Full source-track checks are performed by
lookup/adoption; the completion audit does not recursively revalidate all historical resources.

The core D7 report remains separate from the extended evidence ledgers; this target does not
add prediction scores to the core ranking or change hypothesis weights. Scientific review
and clinical suitability are outside this software-completion criterion.

## Validation status
The configured local case completed the core and extended targets with the final integrity
and identity audit passing. This is not a clean-room resource-acquisition test or proof of
numerical identity to an earlier run. Public tests use synthetic coordinates/identifiers.
