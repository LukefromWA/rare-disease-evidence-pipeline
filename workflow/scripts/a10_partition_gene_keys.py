#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from datetime import datetime, timezone

p = argparse.ArgumentParser()
p.add_argument("--input", required=True)
p.add_argument("--vcf", required=True)
p.add_argument("--named", required=True)
p.add_argument("--unresolved", required=True)
p.add_argument("--audit", required=True)
a = p.parse_args()

source = Path(a.input)
vcf = Path(a.vcf)
named_path = Path(a.named)
unresolved_path = Path(a.unresolved)
audit_path = Path(a.audit)

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def missing(value):
    return value is None or value.strip() in {"", ".", "NA", "NaN"}

inputs = {
    "integrated": source,
    "candidate_vcf": vcf,
    "code": Path(__file__).resolve(),
}
before = {name: sha(path) for name, path in inputs.items()}

with gzip.open(source, "rt", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    columns = reader.fieldnames
    required = {"candidate_id", "gene_symbol", "Feature"}
    if not columns or not required.issubset(columns):
        raise SystemExit("Missing required integrated-evidence columns")
    rows = list(reader)

ids = [row["candidate_id"] for row in rows]
if any(missing(value) for value in ids) or len(set(ids)) != len(ids):
    raise SystemExit("Integrated candidate IDs must be nonblank and unique")

named = [row for row in rows if not missing(row["gene_symbol"])]
unresolved = [row for row in rows if missing(row["gene_symbol"])]

if not named:
    raise SystemExit("No symbol-keyed candidates available for this ranking route")
if len(named) + len(unresolved) != len(rows):
    raise SystemExit("Partition row accounting failed")

targets = {row["candidate_id"]: row["Feature"] for row in unresolved}
annotations_by_id = {}
schema = None

if targets:
    with gzip.open(vcf, "rt") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ,"):
                match = re.search(r'Format:\s*([^"]+)', line)
                if not match:
                    raise SystemExit("Cannot parse original CSQ schema")
                schema = match.group(1).split("|")
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            chrom = re.sub(r"^chr", "", fields[0], flags=re.IGNORECASE)
            candidate = ":".join((chrom, fields[1], fields[3], fields[4]))
            if candidate not in targets:
                continue
            if candidate in annotations_by_id:
                raise SystemExit("Duplicate source VCF record for " + candidate)
            if schema is None:
                raise SystemExit("Missing CSQ schema")

            raw_csq = next(
                (item[4:] for item in fields[7].split(";")
                 if item.startswith("CSQ=")),
                "",
            )
            annotations = [
                dict(zip(schema, item.split("|")))
                for item in raw_csq.split(",") if item
            ]
            selected = [
                item for item in annotations
                if item.get("Feature", "") == targets[candidate]
            ]

            annotations_by_id[candidate] = {
                "selected": selected,
                "all": annotations,
                "raw_csq": raw_csq,
            }

    absent = set(targets) - set(annotations_by_id)
    if absent:
        raise SystemExit("Unresolved candidates absent from VCF: " + repr(sorted(absent)))

extra_columns = [
    "annotation_selected_gene_id",
    "annotation_selected_feature_type",
    "annotation_selected_biotype",
    "annotation_selected_match_count",
    "annotation_selected_records_json",
    "annotation_all_gene_ids_json",
    "annotation_original_csq",
    "gene_key_status",
    "primary_model_status",
]

if set(extra_columns).intersection(columns):
    raise SystemExit("Annotation audit columns already exist in source")

evidence = []
for original in unresolved:
    record = annotations_by_id[original["candidate_id"]]
    selected = record["selected"]
    one = selected[0] if len(selected) == 1 else {}
    gene = one.get("Gene", "")
    all_genes = sorted({
        item.get("Gene", "") for item in record["all"]
        if item.get("Gene", "")
    })

    row = dict(original)
    row.update({
        "annotation_selected_gene_id": gene,
        "annotation_selected_feature_type": one.get("Feature_type", ""),
        "annotation_selected_biotype": one.get("BIOTYPE", ""),
        "annotation_selected_match_count": str(len(selected)),
        "annotation_selected_records_json": json.dumps(selected, sort_keys=True),
        "annotation_all_gene_ids_json": json.dumps(all_genes),
        "annotation_original_csq": record["raw_csq"],
        "gene_key_status": (
            "selected_annotation_ambiguous_or_missing"
            if len(selected) != 1 else
            "ensembl_gene_identified_symbol_missing"
            if gene else
            "no_selected_gene_identifier"
        ),
        "primary_model_status": "not_scored_by_symbol_keyed_primary_route",
    })
    evidence.append(row)

after = {name: sha(path) for name, path in inputs.items()}
if before != after:
    raise SystemExit("Input or code changed during partition")

for path in (named_path, unresolved_path, audit_path):
    path.parent.mkdir(parents=True, exist_ok=True)

def write_table(path, data, fieldnames):
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(data)

with tempfile.TemporaryDirectory(
    prefix=".gene_key_stage_", dir=named_path.parent
) as temporary:
    stage = Path(temporary)
    named_stage = stage / "named.tsv.gz"
    unresolved_stage = stage / "unresolved.tsv.gz"
    audit_stage = stage / "audit.json"

    write_table(named_stage, named, columns)
    write_table(unresolved_stage, evidence, columns + extra_columns)

    manifest = {
        "version": "symbol_key_partition_v1",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "counts": {
            "input_variant_rows": len(rows),
            "symbol_keyed_variant_rows": len(named),
            "symbol_keyed_distinct_genes": len({r["gene_symbol"] for r in named}),
            "retained_symbol_unresolved_variant_rows": len(evidence),
        },
        "inputs": {
            name: {"path": str(path.resolve()), "sha256": before[name]}
            for name, path in inputs.items()
        },
        "outputs": {
            "named_integrated_evidence": {
                "path": str(named_path.resolve()),
                "sha256": sha(named_stage),
            },
            "symbol_unresolved_variant_evidence": {
                "path": str(unresolved_path.resolve()),
                "sha256": sha(unresolved_stage),
            },
        },
        "policy": {
            "blank_symbols_not_a_shared_gene": True,
            "symbols_invented": False,
            "original_integrated_values_changed": False,
            "posterior_values_renormalized": False,
            "unresolved_variants_retained": True,
            "interpretation": (
                "Symbol-unresolved rows remain separate variant evidence. "
                "They are not assigned frozen-model hypothesis weights or "
                "drug scenarios by the symbol-keyed primary route."
            ),
        },
    }
    audit_stage.write_text(json.dumps(manifest, indent=2) + "\n")

    os.replace(named_stage, named_path)
    os.replace(unresolved_stage, unresolved_path)
    os.replace(audit_stage, audit_path)

print(json.dumps(manifest["counts"], indent=2))
