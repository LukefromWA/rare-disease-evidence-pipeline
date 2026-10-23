#!/usr/bin/env bash
# Archive code, configs, logs, results, manifests, environment, and checksums; verify the archive extracts and checksums match.
set -euo pipefail
OUT="${1:-submission_$(date -u +%Y%m%dT%H%M%SZ)}"; mkdir -p "$OUT"
for d in config workflow tests reports/generated reports/output results docs; do [[ -e "$d" ]] && cp -r "$d" "$OUT/" ; done
rm -f "$OUT/config/config.private.yaml"                                   # never ship private absolute paths
{ Rscript -e 'sessionInfo()'; snakemake --version; bcftools --version | head -1; } > "$OUT/ENVIRONMENT.txt" 2>&1 || true
( cd "$OUT" && find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt )
tar -czf "$OUT.tar.gz" "$OUT"; mkdir -p "$OUT.verify"; tar -xzf "$OUT.tar.gz" -C "$OUT.verify"
( cd "$OUT.verify/$OUT" && sha256sum -c SHA256SUMS.txt --quiet ) && echo "freeze OK: $OUT.tar.gz verified"
