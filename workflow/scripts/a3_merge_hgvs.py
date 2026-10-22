#!/usr/bin/env python3
"""A3: exact allele + transcript HGVS merge. Streams the A0 handoff VCF; never matches by gene symbol.
Usage: a3_merge_hgvs.py --vcf V --shortlist S.tsv.gz --out OUT.tsv.gz --audit AUDIT.tsv [--hgvsc HGVSc --hgvsp HGVSp]"""
import argparse, csv, gzip, sys, urllib.parse
from collections import Counter

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--vcf", required=True); ap.add_argument("--shortlist", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--audit", required=True); ap.add_argument("--hgvsc", default="HGVSc"); ap.add_argument("--hgvsp", default="HGVSp")
    a = ap.parse_args()
    with gzip.open(a.shortlist, "rt", newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t"); rows = list(rd); cols = list(rd.fieldnames)
    if "candidate_id" not in cols or "Feature" not in cols:
        sys.exit("A3 STOP: shortlist needs candidate_id and Feature (transcript) columns")
    want = {r["candidate_id"]: r for r in rows}; schema = None; found = {}
    with gzip.open(a.vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                schema = line.split("Format: ")[1].split('"')[0].split("|"); continue
            if line.startswith("#"): continue
            if schema is None: sys.exit("A3 STOP: CSQ header not found before records")
            f = line.rstrip("\n").split("\t", 8)
            cid = f"{f[0].removeprefix('chr')}:{f[1]}:{f[3]}:{f[4]}"
            if cid not in want: continue
            csq = next((x[4:] for x in f[7].split(";") if x.startswith("CSQ=")), None)
            if csq is None: continue
            for tr in csq.split(","):
                vals = tr.split("|"); rec = dict(zip(schema, vals))
                if rec.get("Feature") == want[cid]["Feature"]:
                    found[cid] = (rec.get(a.hgvsc, ""), rec.get(a.hgvsp, "")); break
    def clean_p(x): return urllib.parse.unquote(x.split(":", 1)[1] if ":" in x else x) if x else ""
    def clean_c(x): return x.split(":", 1)[1] if ":" in x else x
    status = Counter()
    for r in rows:
        cid = r["candidate_id"]
        if cid in found:
            hc, hp = found[cid]; r[a.hgvsc] = clean_c(hc); r[a.hgvsp] = clean_p(hp); r["hgvs_merge_status"] = "matched" if (hc or hp) else "matched_no_hgvs"
        else:
            r.setdefault(a.hgvsc, ""); r.setdefault(a.hgvsp, ""); r["hgvs_merge_status"] = "unmatched_transcript_or_no_csq"
        status[r["hgvs_merge_status"]] += 1
    out_cols = cols + [c for c in (a.hgvsc, a.hgvsp, "hgvs_merge_status") if c not in cols]
    with gzip.open(a.out, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    with open(a.audit, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["status", "n"]); [w.writerow([k, v]) for k, v in sorted(status.items())]
        w.writerow(["with_hgvsp", sum(1 for r in rows if r.get(a.hgvsp))]); w.writerow(["n_shortlist", len(rows)])
    if len(rows) != len(want): sys.exit("A3 STOP: duplicate candidate_id in shortlist")
    print(f"A3: {status['matched']} matched, {status['unmatched_transcript_or_no_csq']} unmatched of {len(rows)}", file=sys.stderr)

if __name__ == "__main__": main()
