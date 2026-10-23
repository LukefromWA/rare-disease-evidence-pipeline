#!/usr/bin/env bash
# A0 (existing_vep): validate and stage an already-produced VEP biallelic VCF with full provenance.
# Usage: a0_existing_vep.sh --config C --sample S --vcf V --tbi T --metadata M --outdir D [--raw_vcf R --raw_tbi RT]
set -euo pipefail
while [[ $# -gt 0 ]]; do case "$1" in --config) CONFIG="$2";; --sample) SAMPLE="$2";; --vcf) VCF="$2";; --tbi) TBI="$2";; --metadata) META="$2";; --outdir) OUT="$2";; --raw_vcf) RAW="$2";; --raw_tbi) RAWTBI="$2";; *) echo "unknown arg $1" >&2; exit 2;; esac; shift 2; done
for v in CONFIG SAMPLE VCF TBI META OUT; do [[ -n "${!v:-}" ]] || { echo "missing --${v,,}" >&2; exit 2; }; done
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fail() { echo "A0 STOP: $*" >&2; exit 1; }
[[ -s "$VCF" ]] || fail "VCF absent or empty: $VCF"
[[ -s "$TBI" ]] || fail "VCF index absent: $TBI"
[[ -s "$META" ]] || fail "VEP run metadata absent (existing_vep needs provenance): $META"
gzip -t "$VCF" || fail "VCF is not valid gzip: $VCF"
HEADER="$(gzip -cd "$VCF" | grep '^#' | head -n 100000)"
grep -q '^##INFO=<ID=CSQ' <<<"$HEADER" || fail "no CSQ INFO header (not VEP-annotated)"
CHROMLINE="$(grep '^#CHROM' <<<"$HEADER")"
grep -qw "$SAMPLE" <<<"$CHROMLINE" || echo "A0 note: pipeline sample label is not a VCF column name; A1 validates vcf.sample_id or vcf.sample_index." >&2
mkdir -p "$OUT"
cp -f "$VCF" "$OUT/annotated.biallelic.vcf.gz"; cp -f "$TBI" "$OUT/annotated.biallelic.vcf.gz.tbi"
cmp -s "$VCF" "$OUT/annotated.biallelic.vcf.gz" || fail "staged copy differs from source"
Rscript "$HERE/a0_manifest.R" --config "$CONFIG" --mode existing_vep --sample "$SAMPLE" --vcf "$OUT/annotated.biallelic.vcf.gz" \
  --tbi "$OUT/annotated.biallelic.vcf.gz.tbi" --source_vcf "$VCF" --metadata "$META" --outdir "$OUT" --raw_vcf "${RAW:-NA}" --raw_tbi "${RAWTBI:-NA}"
