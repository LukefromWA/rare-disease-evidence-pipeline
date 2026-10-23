#!/usr/bin/env python3
"""Layer 3: targeted openFDA label lookup for the drugs that D2 actually surfaced.

Reads a drug table (D2 output), queries only distinct drugs, caches every terminal API response,
and writes label evidence, a human review worksheet and a D5-compatible pediatric_safety_map.
It never infers safety: statuses stay 'unknown' unless a human-reviewed override cites a label ID.
Only public drug names are sent; no patient or candidate data. Standard library + system curl.
"""
import argparse, csv, gzip, hashlib, json, os, re, subprocess, sys, tempfile, time
from pathlib import Path
from urllib.parse import quote
import datetime as dt

VERSION = "track2_openfda_targeted_labels_v1"
API = "https://api.fda.gov/drug/label.json"
SECTIONS = ["pediatric_use", "use_in_specific_populations", "contraindications", "boxed_warning",
            "warnings_and_cautions", "warnings", "indications_and_usage"]
PED = {"supportive", "limited", "unknown", "contraindicated"}
CONTRA = {"none", "relative", "absolute", "unknown"}
PAGE = 100

class TransportError(RuntimeError):
    pass

def sha_bytes(b): return hashlib.sha256(b).hexdigest()
def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def read_tsv(p):
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt", encoding="utf-8-sig", newline="") as f:
        r = csv.DictReader(f, delimiter="\t"); cols = r.fieldnames or []
        rows = list(r)
    if not cols or len(cols) != len(set(cols)): raise ValueError("Missing/duplicate header: " + str(p))
    if any(None in x for x in rows): raise ValueError("Ragged TSV row: " + str(p))
    return cols, rows

def write_tsv(p, rows, cols):
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "wt", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="raise")
        w.writeheader(); w.writerows(rows)
    if len(read_tsv(p)[1]) != len(rows): raise ValueError("Read-back failed: " + str(p))

def search_url(search, limit, skip):
    return "%s?search=%s&limit=%d&skip=%d" % (API, quote(search, safe=":"), limit, skip)

def curl_transport(curl, key):
    def get(url):
        full = url + ("&api_key=" + key if key else "")
        with tempfile.NamedTemporaryFile(delete=False) as t: body = t.name
        try:
            cmd = [curl, "-K", "-", "--proto", "=https", "--proto-redir", "=https", "--http1.1", "--silent", "--show-error",
                   "--location", "--connect-timeout", "30", "--max-time", "120", "--output", body, "--write-out", "%{http_code}"]
            # URL (and key) go through stdin config so the key never appears in the process list.
            r = subprocess.run(cmd, input='url = "%s"\n' % full, text=True, capture_output=True)
            code = int((r.stdout or "0").strip() or 0)
            if r.returncode != 0 and code == 0: raise TransportError("curl exit %d" % r.returncode)
            return code, Path(body).read_bytes()
        finally:
            Path(body).unlink(missing_ok=True)
    return get

class Client:
    def __init__(self, cache, transport, delay, sleep=time.sleep):
        self.cache, self.transport, self.delay, self.sleep = Path(cache), transport, delay, sleep
        self.cache.mkdir(parents=True, exist_ok=True); self.hits = 0; self.calls = 0
        self.last_updated = set(); self.cache_keys = []

    def get(self, url):
        key = sha_bytes(url.encode()); body_p = self.cache / (key + ".json"); meta_p = self.cache / (key + ".meta.json")
        self.cache_keys.append(key)
        if body_p.exists() and meta_p.exists():
            meta = json.loads(meta_p.read_text()); body = body_p.read_bytes()
            if meta["url"] != url or meta["body_sha256"] != sha_bytes(body): raise ValueError("Corrupt cache entry: " + key)
            self.hits += 1
            return meta["http_code"], json.loads(body)
        last = None
        for attempt in range(4):
            try:
                code, body = self.transport(url); self.calls += 1
            except TransportError as e:
                last = str(e); code, body = 0, b""
            if code == 200 or (code == 404 and b"NOT_FOUND" in body):
                obj = json.loads(body)
                if code == 200 and (not isinstance(obj.get("results"), list) or "meta" not in obj):
                    raise ValueError("Unexpected openFDA response schema")
                body_p.write_bytes(body)
                meta_p.write_text(json.dumps({"url": url, "http_code": code, "body_sha256": sha_bytes(body),
                                              "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat()}, indent=2))
                if code == 200 and obj["meta"].get("last_updated"): self.last_updated.add(obj["meta"]["last_updated"])
                self.sleep(self.delay)
                return code, obj
            last = last or "HTTP %s" % code
            if code not in (0, 429, 500, 502, 503, 504): break
            self.sleep(2 ** attempt)
        raise TransportError(last or "request failed")

def quote_phrase(s): return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'

def fetch_labels(client, search, cap):
    labels, total, skip = [], None, 0
    while True:
        limit = min(PAGE, cap - len(labels))
        code, obj = client.get(search_url(search, limit, skip))
        if code == 404: return [], 0
        res = obj["meta"].get("results", {}); total = int(res.get("total", len(obj["results"])))
        labels += obj["results"]; skip += len(obj["results"])
        if not obj["results"] or len(labels) >= min(total, cap): return labels, total

def lookup_drug(client, name, cap):
    exact = 'openfda.generic_name.exact:' + quote_phrase(name.upper())
    labels, total = fetch_labels(client, exact, cap)
    if labels: return "matched_exact", labels, total
    labels, total = fetch_labels(client, "openfda.generic_name:" + quote_phrase(name), cap)
    if labels: return "matched_phrase_review", labels, total
    return "no_match", [], 0

def text_of(v): return "\n".join(v) if isinstance(v, list) else str(v)
def pipe(v): return "|".join(v) if isinstance(v, list) else (v or "")

def load_overrides(path, evidence, drugs):
    if not path: return {}
    cols, rows = read_tsv(path)
    need = ["drug_id", "pediatric_evidence_status", "contraindication_status", "route_feasible", "evidence_label_id", "reviewer_note"]
    if any(c not in cols for c in need): raise ValueError("Overrides need columns: " + ", ".join(need))
    out = {}
    for r in rows:
        d = r["drug_id"].strip()
        if not d: continue
        if d in out: raise ValueError("Duplicate override for " + d)
        if d not in drugs: raise ValueError("Override for drug not in input: " + d)
        ped, con, route = r["pediatric_evidence_status"].strip().lower(), r["contraindication_status"].strip().lower(), r["route_feasible"].strip().upper()
        if ped not in PED or con not in CONTRA or route not in ("TRUE", "FALSE"): raise ValueError("Out-of-vocabulary override for " + d)
        if not r["reviewer_note"].strip(): raise ValueError("Reviewer note required for " + d)
        if (d, r["evidence_label_id"].strip()) not in evidence: raise ValueError("evidence_label_id must be a label retrieved for %s" % d)
        out[d] = {**r, "pediatric_evidence_status": ped, "contraindication_status": con, "route_feasible": route}
    return out

def main(argv=None, transport=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--drugs", required=True); p.add_argument("--out-dir", required=True); p.add_argument("--cache-dir", required=True)
    p.add_argument("--drug-id-column", default="drug_id"); p.add_argument("--drug-name-column", default="drug_name")
    p.add_argument("--review-overrides"); p.add_argument("--max-drugs", type=int, default=200)
    p.add_argument("--max-labels-per-drug", type=int, default=100); p.add_argument("--delay", type=float, default=0.6)
    p.add_argument("--curl", default="/usr/bin/curl")
    a = p.parse_args(argv)
    out = Path(a.out_dir); names = ["drug_lookup_status.tsv", "label_evidence.tsv.gz", "label_sections.tsv.gz", "review_worksheet.tsv",
                                    "pediatric_safety_map.tsv.gz", "run_manifest.json", "checksums.sha256"]
    out.mkdir(parents=True, exist_ok=True)
    if any((out / n).exists() for n in names): raise ValueError("Refusing to overwrite existing outputs in " + str(out))
    cols, rows = read_tsv(a.drugs)
    if a.drug_id_column not in cols or a.drug_name_column not in cols: raise ValueError("Drug table missing id/name columns")
    drugs = {}
    for r in rows:
        d = r[a.drug_id_column].strip()
        if not d: continue
        n = r[a.drug_name_column].strip()
        if d in drugs and drugs[d] and n and drugs[d] != n: raise ValueError("Conflicting names for " + d)
        drugs[d] = drugs.get(d) or n
    if not drugs: raise ValueError("No drugs supplied")
    if len(drugs) > a.max_drugs: raise ValueError("%d distinct drugs exceeds --max-drugs %d; narrow selection or raise it deliberately" % (len(drugs), a.max_drugs))
    key = os.environ.get("OPENFDA_API_KEY")
    if not key: print("WARNING: OPENFDA_API_KEY not set; requests may be rate-limited or refused.", file=sys.stderr)
    client = Client(a.cache_dir, transport or curl_transport(a.curl, key), a.delay)
    status, ev, sec, evidence_ids = [], [], [], set()
    for d, name in sorted(drugs.items()):
        if not name:
            status.append({"drug_id": d, "drug_name": "", "status": "missing_name", "labels_fetched": 0, "labels_total": 0, "result_complete": "NA"}); continue
        try:
            st, labels, total = lookup_drug(client, name, a.max_labels_per_drug)
        except TransportError as e:
            status.append({"drug_id": d, "drug_name": name, "status": "request_failed", "labels_fetched": 0, "labels_total": "NA", "result_complete": "NA", "note": str(e)}); continue
        status.append({"drug_id": d, "drug_name": name, "status": st, "labels_fetched": len(labels), "labels_total": total,
                       "result_complete": "yes" if len(labels) >= total else "no", "note": ""})
        for L in labels:
            lid = L.get("id", ""); o = L.get("openfda") or {}
            if not lid: raise ValueError("Label without id for " + d)
            evidence_ids.add((d, lid))
            row = {"drug_id": d, "label_id": lid, "set_id": L.get("set_id", ""), "version": L.get("version", ""),
                   "effective_time": L.get("effective_time", ""), "match_type": st, "product_type": pipe(o.get("product_type")),
                   "route": pipe(o.get("route")), "generic_name": pipe(o.get("generic_name")), "brand_name": pipe(o.get("brand_name")),
                   "unii": pipe(o.get("unii")), "application_number": pipe(o.get("application_number"))}
            for s in SECTIONS:
                row["has_" + s] = "yes" if L.get(s) else "no"
                if L.get(s): sec.append({"drug_id": d, "label_id": lid, "section": s, "text": text_of(L[s])})
            ev.append(row)
    overrides = load_overrides(a.review_overrides, evidence_ids, drugs)
    release = ";".join(sorted(client.last_updated)) or "unknown"
    safety, work = [], []
    for d, name in sorted(drugs.items()):
        o = overrides.get(d); mine = [e for e in ev if e["drug_id"] == d]
        safety.append({"drug_id": d, "pediatric_evidence_status": o["pediatric_evidence_status"] if o else "unknown",
            "contraindication_status": o["contraindication_status"] if o else "unknown", "route_feasible": o["route_feasible"] if o else "",
            "review_state": "reviewed" if o else "unreviewed", "evidence_label_id": o["evidence_label_id"] if o else "",
            "reviewer_note": o["reviewer_note"] if o else "", "source_database": "openFDA drug label API (human review where review_state=reviewed)",
            "source_release": release})
        ex = next((e["label_id"] for e in mine if e["has_pediatric_use"] == "yes"), mine[0]["label_id"] if mine else "")
        work.append({"drug_id": d, "drug_name": name, "n_labels": len(mine), "any_pediatric_use_section": "yes" if any(e["has_pediatric_use"] == "yes" for e in mine) else "no",
            "any_boxed_warning": "yes" if any(e["has_boxed_warning"] == "yes" for e in mine) else "no",
            "any_contraindications_section": "yes" if any(e["has_contraindications"] == "yes" for e in mine) else "no", "example_label_id": ex,
            "pediatric_evidence_status": "", "contraindication_status": "", "route_feasible": "", "evidence_label_id": "", "reviewer_note": ""})
    ev_cols = ["drug_id", "label_id", "set_id", "version", "effective_time", "match_type", "product_type", "route", "generic_name", "brand_name", "unii", "application_number"] + ["has_" + s for s in SECTIONS]
    write_tsv(out / "drug_lookup_status.tsv", status, ["drug_id", "drug_name", "status", "labels_fetched", "labels_total", "result_complete", "note"])
    write_tsv(out / "label_evidence.tsv.gz", ev, ev_cols)
    write_tsv(out / "label_sections.tsv.gz", sec, ["drug_id", "label_id", "section", "text"])
    write_tsv(out / "review_worksheet.tsv", work, list(work[0].keys()))
    write_tsv(out / "pediatric_safety_map.tsv.gz", safety, list(safety[0].keys()))
    counts = {}
    for s in status: counts[s["status"]] = counts.get(s["status"], 0) + 1
    man = {"adapter_version": VERSION, "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "script_sha256": sha(__file__),
        "drugs_input_sha256": sha(a.drugs), "overrides_sha256": sha(a.review_overrides) if a.review_overrides else None,
        "distinct_drugs": len(drugs), "status_counts": counts, "reviewed_drugs": len(overrides), "api_requests_made": client.calls,
        "cache_hits": client.hits, "api_meta_last_updated": release, "cache_entries_used": sorted(set(client.cache_keys)),
        "limitations": ["Only drug names were sent to openFDA; no patient or candidate data.",
            "matched_phrase_review can include combination products or salts and needs identity review.",
            "Section presence is not approval, pediatric suitability or safety. Unreviewed drugs remain 'unknown' with no route_feasible, which D5 treats as not_configured.",
            "Overrides must be human-reviewed and cite a retrieved label_id. This is not medical advice."]}
    (out / "run_manifest.json").write_text(json.dumps(man, indent=2) + "\n", encoding="utf-8")
    (out / "checksums.sha256").write_text("".join(sha(out / n) + "  " + n + "\n" for n in names[:-1]))
    print("LABEL EVIDENCE COMPLETE:", out); print(json.dumps({k: man[k] for k in ("distinct_drugs", "status_counts", "reviewed_drugs", "api_requests_made", "cache_hits")}, indent=2))

if __name__ == "__main__":
    try: main()
    except (ValueError, OSError, TransportError) as e:
        print("STOP:", e, file=sys.stderr); sys.exit(1)
