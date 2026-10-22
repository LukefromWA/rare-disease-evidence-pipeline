#!/usr/bin/env bash
# A0 (rerun_vep): normalize/split with bcftools, annotate with VEP from pinned assets, write the resolved command.
# Usage: a0_normalize_and_vep.sh --config C --sample S --vcf V --tbi T --fasta F --fai FI --cache_dir D --outdir O [--threads N --cache_version 116]
set -euo pipefail
THREADS=4; CACHE_VERSION=116
while [[ $# -gt 0 ]]; do case "$1" in --config) CONFIG="$2";; --sample) SAMPLE="$2";; --vcf) VCF="$2";; --tbi) TBI="$2";; --fasta) FASTA="$2";; --fai) FAI="$2";; --cache_dir) CACHE="$2";; --outdir) OUT="$2";; --threads) THREADS="$2";; --cache_version) CACHE_VERSION="$2";; *) echo "unknown arg $1" >&2; exit 2;; esac; shift 2; done
fail() { echo "A0 STOP: $*" >&2; exit 1; }
for v in CONFIG SAMPLE VCF TBI FASTA FAI CACHE OUT; do [[ -n "${!v:-}" && "${!v}" != "null" && "${!v}" != "NA" ]] || fail "required asset not configured (${v,,}); rerun_vep fails loudly instead of inventing a path"; done
for f in "$VCF" "$TBI" "$FASTA" "$FAI"; do [[ -s "$f" ]] || fail "missing file: $f"; done; [[ -d "$CACHE" ]] || fail "VEP cache dir absent: $CACHE"
command -v bcftools >/dev/null || fail "bcftools not on PATH"; command -v vep >/dev/null || fail "vep not on PATH"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; mkdir -p "$OUT"
NORM="$OUT/normalized.biallelic.vcf.gz"; ANN="$OUT/annotated.biallelic.vcf.gz"
CMD1=(bcftools norm -f "$FASTA" -m -any --check-ref w "$VCF" -Oz -o "$NORM")
CMD2=(vep --input_file "$NORM" --output_file "$ANN" --vcf --compress_output bgzip --assembly GRCh38 --species homo_sapiens --cache --offline --cache_version "$CACHE_VERSION" --dir_cache "$CACHE" --fasta "$FASTA" --hgvs --symbol --canonical --mane --max_af --af_gnomade --af_gnomadg --fork "$THREADS" --force_overwrite)
{ printf '%q ' "${CMD1[@]}"; echo; printf '%q ' "${CMD2[@]}"; echo; } > "$OUT/resolved_vep_command.txt"      # recorded verbatim
"${CMD1[@]}"; tabix -f -p vcf "$NORM"; "${CMD2[@]}"; tabix -f -p vcf "$ANN"
Rscript "$HERE/a0_manifest.R" --config "$CONFIG" --mode rerun_vep --sample "$SAMPLE" --vcf "$ANN" --tbi "$ANN.tbi" --source_vcf "$VCF" --metadata "$OUT/resolved_vep_command.txt" \
  --outdir "$OUT" --raw_vcf "$VCF" --raw_tbi "$TBI" --fasta "$FASTA" --fai "$FAI"
