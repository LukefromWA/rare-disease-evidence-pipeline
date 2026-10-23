#!/usr/bin/env python3
"""Layer 2 bridge: Track 1 layer-A gene-level candidates -> D0 input (C4-format hypothesis set).

This is a DECLARED A-only amendment. qj here is a normalized layer-A posterior over the supplied rows,
not a frozen B/C-model probability. Standard library only; no network; no case names hardcoded.
Output columns follow the build-sheet C4 contract: gene, variant, qj, mechanism, pm (+ audit columns).
"""
import argparse, csv, gzip, hashlib, json, math, re, sys
from pathlib import Path
import datetime as dt

VERSION = "track2_a10_to_d0_bridge_v1"
LOF = re.compile(r"stop_gained|frameshift|splice_donor|splice_acceptor|start_lost|stop_lost")
MIS = re.compile(r"missense|inframe")
DEFAULT_MECHANISMS = ["loss_of_function", "gain_of_function", "dominant_negative", "unknown"]

def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def read_tsv(p):
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt", encoding="utf-8-sig", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        cols = r.fieldnames or []
        if not cols or len(cols) != len(set(cols)):
            raise ValueError("Missing/duplicate header in " + str(p))
        rows = list(r)
    if any(None in x for x in rows):
        raise ValueError("Ragged TSV row in " + str(p))
    return cols, rows

def consequence_class(consequence, impact):
    """Mirrors the build-sheet C4 chunk as extracted; verify against your hand-coded C4."""
    c = (consequence or "").strip().lower()
    i = (impact or "").strip().upper()
    if c:
        if LOF.search(c): return "lof_like"
        if MIS.search(c): return "missense"
        return "other"
    return "lof_like" if i == "HIGH" else "missense" if i == "MODERATE" else "other"

def load_prior(path, classes_present):
    cols, rows = read_tsv(path)
    need = {"consequence_class", "mechanism", "pm"}
    if not need <= set(cols):
        raise ValueError("Prior table needs columns: consequence_class, mechanism, pm")
    table = {}
    for r in rows:
        try:
            pm = float(r["pm"])
        except ValueError:
            raise ValueError("Non-numeric pm for %s/%s: fill the declared prior (decision D-3) or use --prior-mode uniform_neutral" % (r["consequence_class"], r["mechanism"]))
        if not math.isfinite(pm) or pm < 0 or not r["mechanism"].strip():
            raise ValueError("Invalid mechanism/pm in prior table")
        key = r["consequence_class"].strip()
        if r["mechanism"].strip() in dict(table.get(key, [])):
            raise ValueError("Duplicate prior row: %s/%s" % (key, r["mechanism"]))
        table.setdefault(key, []).append((r["mechanism"].strip(), pm))
    for k, v in table.items():
        if abs(sum(p for _, p in v) - 1) > 1e-8:
            raise ValueError("Prior rows for class %r must sum to 1" % k)
    out = {}
    for c in classes_present:
        if c in table: out[c] = table[c]
        elif "default" in table: out[c] = table["default"]
        else: raise ValueError("No prior for consequence class %r (add it or a 'default' class)" % c)
    return out

def select(rows, k, ties):
    ordered = sorted(rows, key=lambda r: (-r["_post"], r["_gene"], r["_variant"]))
    if k is None or len(ordered) <= k:
        return ordered
    cutoff = ordered[k - 1]["_post"]
    return [r for r in ordered if r["_post"] >= cutoff] if ties == "include" else ordered[:k]

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--a10", required=True)
    p.add_argument("--evidence-matrix", help="Optional TSV with candidate_id + consequence for consequence classes")
    p.add_argument("--out", required=True)
    p.add_argument("--gene-column", default="gene_symbol")
    p.add_argument("--variant-column", default="best_candidate_id")
    p.add_argument("--impact-column", default="best_IMPACT")
    p.add_argument("--posterior-column", default="posterior_informed")
    p.add_argument("--evidence-id-column", default="candidate_id")
    p.add_argument("--evidence-consequence-column", default="consequence")
    p.add_argument("--selection", choices=["full", "top_k"], default="full")
    p.add_argument("--top-k", type=int)
    p.add_argument("--tie-policy", choices=["include", "truncate"], default="include")
    p.add_argument("--prior-mode", choices=["declared", "uniform_neutral"], required=True)
    p.add_argument("--mechanism-prior")
    p.add_argument("--mechanisms", default=",".join(DEFAULT_MECHANISMS))
    a = p.parse_args(argv)
    out = Path(a.out)
    manifest_path = out.with_name(out.name.replace(".tsv.gz", "") + ".manifest.json")
    if out.exists() or manifest_path.exists():
        raise ValueError("Refusing to overwrite existing output")
    if a.selection == "top_k" and (a.top_k is None or a.top_k < 1):
        raise ValueError("--selection top_k needs --top-k >= 1")
    if a.prior_mode == "declared" and not a.mechanism_prior:
        raise ValueError("--prior-mode declared needs --mechanism-prior")
    cols, rows = read_tsv(a.a10)
    need = [a.gene_column, a.variant_column, a.posterior_column]
    if any(c not in cols for c in need):
        raise ValueError("A10 missing columns: " + ", ".join(c for c in need if c not in cols))
    seen = set(); recs = []
    for r in rows:
        gene, var = r[a.gene_column].strip(), r[a.variant_column].strip()
        if not gene or not var:
            raise ValueError("Blank gene or variant identifier in A10")
        if var in seen:
            raise ValueError("Duplicate variant identifier: " + var)
        seen.add(var)
        try:
            post = float(r[a.posterior_column])
        except ValueError:
            raise ValueError("Non-numeric posterior for " + var)
        if not math.isfinite(post) or post < 0:
            raise ValueError("Posterior must be finite and non-negative: " + var)
        recs.append({**r, "_gene": gene, "_variant": var, "_post": post})
    if not recs: raise ValueError("Empty A10 table")
    total = sum(r["_post"] for r in recs)
    if total <= 0: raise ValueError("A10 posterior mass is zero")
    kept = select(recs, a.top_k if a.selection == "top_k" else None, a.tie_policy)
    kept_mass = sum(r["_post"] for r in kept)
    cons = {}
    if a.evidence_matrix:
        ec, er = read_tsv(a.evidence_matrix)
        if a.evidence_id_column not in ec or a.evidence_consequence_column not in ec:
            raise ValueError("Evidence matrix missing id/consequence columns")
        cons = {x[a.evidence_id_column]: x[a.evidence_consequence_column] for x in er}
    for r in kept:
        r["_class"] = consequence_class(cons.get(r["_variant"]), r.get(a.impact_column))
    present = sorted({r["_class"] for r in kept})
    mechs = [m.strip() for m in a.mechanisms.split(",") if m.strip()]
    if a.prior_mode == "declared":
        prior = load_prior(a.mechanism_prior, present)
    else:
        print("WARNING: uniform_neutral prior is NOT a declared biological prior (decision D-3). Label it so in any report.", file=sys.stderr)
        prior = {c: [(m, 1.0 / len(mechs)) for m in mechs] for c in present}
    lines = []
    for r in kept:
        qj = r["_post"] / kept_mass if kept_mass > 0 else 1.0 / len(kept)
        for mech, pm in prior[r["_class"]]:
            lines.append([r["_gene"], r["_variant"], "%.12g" % qj, mech, "%.12g" % pm, r["_class"], a.posterior_column, a.prior_mode])
    uq = {}
    for g, v, q, *_ in lines: uq[v] = float(q)
    if abs(sum(uq.values()) - 1) > 1e-8: raise ValueError("qj does not sum to 1 over selected hypotheses")
    pmsum = {}
    for g, v, q, m, pm, *_ in lines: pmsum[v] = pmsum.get(v, 0) + float(pm)
    if any(abs(x - 1) > 1e-8 for x in pmsum.values()): raise ValueError("pm does not sum to 1 within a variant")
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "variant", "qj", "mechanism", "pm", "consequence_class", "posterior_source", "prior_mode"])
        w.writerows(lines)
    with gzip.open(out, "rt", encoding="utf-8") as f:
        if sum(1 for _ in f) - 1 != len(lines): raise ValueError("Read-back row count failed")
    manifest = {"adapter_version": VERSION, "built_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "output_sha256": sha(out), "a10_sha256": sha(a.a10), "evidence_matrix_sha256": sha(a.evidence_matrix) if a.evidence_matrix else None,
        "mechanism_prior_sha256": sha(a.mechanism_prior) if a.mechanism_prior and a.prior_mode == "declared" else None,
        "script_sha256": sha(__file__), "posterior_column": a.posterior_column,
        "selection": {"mode": a.selection, "top_k": a.top_k, "tie_policy": a.tie_policy, "input_rows": len(recs),
                      "selected_hypotheses": len(kept), "truncated_mass_fraction": round(1 - kept_mass / total, 12)},
        "prior": {"mode": a.prior_mode, "is_declared_biological_prior": a.prior_mode == "declared", "mechanisms": mechs if a.prior_mode != "declared" else None},
        "counts": {"hypotheses": len(uq), "scenario_rows": len(lines)},
        "limitations": ["A-only amendment: qj is a normalized layer-A posterior, not a frozen B/C model probability.",
            "With top_k, qj is renormalized over the kept rows; truncated_mass_fraction reports what was dropped.",
            "Consequence-class rules mirror the build-sheet C4 chunk as extracted; verify against the hand-coded C4.",
            "Mechanism names must match those used by the D1/D3 code. No clinical conclusions."]}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print("D0 INPUT WRITTEN:", out)
    print(json.dumps({"hypotheses": len(uq), "scenario_rows": len(lines), "selection": manifest["selection"], "prior_mode": a.prior_mode}, indent=2))

if __name__ == "__main__":
    try: main()
    except (ValueError, OSError) as e:
        print("STOP:", e, file=sys.stderr); sys.exit(1)
