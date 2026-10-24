#!/usr/bin/env bash
set -Eeuo pipefail
POLICY=legacy_any_csq_token_position_v1
INPUT= OUTDIR= PROGRESS=100000
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT="${2:?}"; shift 2;;
    --outdir) OUTDIR="${2:?}"; shift 2;;
    --policy) POLICY="${2:?}"; shift 2;;
    --progress-every) PROGRESS="${2:?}"; shift 2;;
    *) echo "Unknown argument: $1" >&2; exit 2;;
  esac
done
[[ -f "$INPUT" && -n "$OUTDIR" ]] || { echo 'Required: --input VCF --outdir DIR' >&2; exit 2; }
[[ "$POLICY" == legacy_any_csq_token_position_v1 ]] || { echo 'Unsupported candidate policy' >&2; exit 2; }
[[ "$PROGRESS" =~ ^[1-9][0-9]*$ ]] || { echo 'progress-every must be a positive integer' >&2; exit 2; }
for tool in bcftools awk python3 sha256sum gzip; do
  command -v "$tool" >/dev/null || { echo "Missing command: $tool" >&2; exit 2; }
done
INPUT=$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve())' "$INPUT")
mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd -P)
for product in candidates.vcf.gz candidates.vcf.gz.tbi candidate_manifest.json; do
  [[ ! -e "$OUTDIR/$product" ]] || { echo "Refusing existing product: $OUTDIR/$product" >&2; exit 1; }
done
STAGE=$(mktemp -d "$OUTDIR/.candidate_stage.XXXXXXXX")
trap 'rm -rf -- "$STAGE"' EXIT
sha256sum "$INPUT" > "$STAGE/before.sha256"
printf 'Candidate policy: %s\n' "$POLICY" >&2
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ\n' "$INPUT" \
| awk -F'\t' -v counts="$STAGE/counts.tsv" -v progress="$PROGRESS" '
  BEGIN { OFS="\t" }
  {
    if (toupper($5) ~ /(^|[|,])(HIGH|MODERATE)([|,]|$)/) {
      print $1, $2
      kept++
    }
    if (NR % progress == 0) {
      printf "Selection progress: scanned %d records; selected %d records\n", NR, kept > "/dev/stderr"
      fflush("/dev/stderr")
    }
  }
  END {
    print NR+0, kept+0 > counts
    close(counts)
    printf "Selection complete: scanned %d records; selected %d records\n", NR, kept > "/dev/stderr"
  }
' > "$STAGE/positions.tsv"
[[ -s "$STAGE/positions.tsv" ]] || { echo 'No candidates selected; no A1 input published' >&2; exit 1; }
bcftools view --no-version --targets-file "$STAGE/positions.tsv" \
  --output-type z --output-file "$STAGE/candidates.vcf.gz" "$INPUT"
bcftools index --tbi "$STAGE/candidates.vcf.gz"
gzip -t "$STAGE/candidates.vcf.gz"
bcftools view --no-header "$STAGE/candidates.vcf.gz" | awk 'END { print NR+0 }' > "$STAGE/subset_count.txt"
python3 - "$INPUT" "$OUTDIR" "$STAGE" "$POLICY" "${BASH_SOURCE[0]}" <<'PYMETA'
import datetime, hashlib, json, pathlib, subprocess, sys
source, outdir, stage = map(pathlib.Path, sys.argv[1:4])
policy = sys.argv[4]
def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()
before = (stage/'before.sha256').read_text().split()[0]
if sha(source) != before:
    raise SystemExit('Source VCF changed during candidate preparation')
records, selected = map(int, (stage/'counts.tsv').read_text().split())
kept = int((stage/'subset_count.txt').read_text())
report = dict(status='complete', created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    policy=policy, selector_implementation='boundary_regex_casefold_v1',
    selection_key=['CHROM','POS'],
    selection_semantics='HIGH or MODERATE token in any CSQ field; retain every record at selected positions',
    input_vcf=str(source), input_vcf_sha256=before,
    subset_vcf=str(outdir/'candidates.vcf.gz'), subset_vcf_sha256=sha(stage/'candidates.vcf.gz'),
    subset_index_sha256=sha(stage/'candidates.vcf.gz.tbi'),
    producer_sha256=sha(pathlib.Path(sys.argv[5])),
    input_records=records, selected_input_records=selected,
    kept_records=kept, records_dropped=records-kept,
    af_filter_applied=False, hard_qc_filters_applied=[], contig_filter_applied=False,
    annotation_retention='full_INFO_FORMAT_CSQ',
    bcftools_version=subprocess.check_output(['bcftools','--version'], text=True).splitlines()[0])
(stage/'candidate_manifest.json').write_text(json.dumps(report, indent=2)+'\n')
PYMETA
mv -- "$STAGE/candidates.vcf.gz" "$OUTDIR/candidates.vcf.gz"
mv -- "$STAGE/candidates.vcf.gz.tbi" "$OUTDIR/candidates.vcf.gz.tbi"
mv -- "$STAGE/candidate_manifest.json" "$OUTDIR/candidate_manifest.json"
printf 'Candidate input published: %s/candidates.vcf.gz\n' "$OUTDIR" >&2
