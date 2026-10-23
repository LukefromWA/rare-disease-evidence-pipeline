#!/usr/bin/env python3
import argparse
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import re
import tempfile

SOURCE_FIELDS = [
    "source_id", "source_type", "title", "authors", "year",
    "publication", "doi", "url", "snapshot", "snapshot_sha256",
]
CLAIM_FIELDS = [
    "claim_id", "source_id", "compound", "evidence_class",
    "source_gene", "population_or_model", "outcome", "claim_text",
    "statistics_json", "source_locator", "source_excerpt",
    "limitations", "dependence_group", "extraction_method", "review_status",
]
LINK_FIELDS = [
    "claim_id", "hypothesis_id", "gene", "variant",
    "relationship", "selection_witness", "q_j_raw",
    "applicability_review_status", "case_limitations",
]
REVIEW_FIELDS = [
    "claim_id", "source_id", "evidence_class",
    "required_review", "status",
]
CLAIM_SOURCE_FIELDS = [
    "claim_id", "source_id", "relationship", "review_status",
]
SNAPSHOT_FIELDS = [
    "source_id", "snapshot", "snapshot_sha256",
]
SOURCE_REVIEW_FIELDS = [
    "source_id", "role", "required_review", "status",
]

def sha(data):
    return hashlib.sha256(data).hexdigest()

def canonical_doi(value):
    value = value.strip().lower()
    for prefix in ["https://doi.org/", "http://doi.org/", "doi:"]:
        if value.startswith(prefix):
            value = value[len(prefix):].strip()
    return value

def combine_sources(base_rows, added_rows):
    sources = {}
    snapshots = {}
    for row in base_rows + added_rows:
        source_id = row["source_id"]
        if not source_id:
            raise ValueError("Blank source ID")
        if source_id in sources:
            previous = sources[source_id]
            first = canonical_doi(previous["doi"])
            second = canonical_doi(row["doi"])
            if first and second:
                if first != second:
                    raise ValueError("Same source ID has conflicting DOIs")
            elif previous["title"].strip() != row["title"].strip():
                raise ValueError("Source identity requires manual reconciliation")
        else:
            sources[source_id] = row.copy()

        key = (
            source_id, row["snapshot"], row["snapshot_sha256"]
        )
        snapshots[key] = dict(zip(SNAPSHOT_FIELDS, key))

    return list(sources.values()), list(snapshots.values())

def combine_claims(base_rows, added_rows):
    claims = {}
    for row in base_rows + added_rows:
        claim_id = row["claim_id"]
        if claim_id in claims and claims[claim_id] != row:
            raise ValueError("Conflicting records share a claim ID")
        claims[claim_id] = row.copy()
    return list(claims.values())

def mechanism_relation(compound):
    if compound.strip().lower() == "g418":
        return "biochemical_comparator_only_not_candidate_nomination"
    return "translation_mechanism_rationale_not_case_functional_validation"

def write_table(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fields,
            delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing table header: {path}")
        rows = list(reader)
    if any(
        None in row or any(value is None for value in row.values())
        for row in rows
    ):
        raise ValueError(f"Malformed table: {path}")
    return rows

def ris_records(value):
    records = {}
    for entry in re.split(r"(?m)(?=^TY  - )", value):
        if not entry.strip():
            continue
        matches = re.findall(r"(?m)^ID  - (.+)$", entry)
        if len(matches) != 1:
            raise ValueError("RIS record has no unique ID")
        source_id = matches[0].strip()
        if source_id in records:
            raise ValueError("Duplicate RIS ID within a reference bundle")
        records[source_id] = entry.strip()
    return records

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-ledger", required=True)
    parser.add_argument("--augmentation", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    base = Path(args.base_ledger)
    augmentation = Path(args.augmentation)
    tracked = {}

    def track(path):
        path = Path(path)
        actual = sha(path.read_bytes())
        old = tracked.get(str(path))
        if old is not None and old != actual:
            raise ValueError("Input changed during reading")
        tracked[str(path)] = actual
        return actual

    def load_json(path):
        track(path)
        return json.loads(Path(path).read_text(encoding="utf-8"))

    def verify_bundle(directory):
        manifest = load_json(directory / "manifest.json")
        hashes = manifest.get("output_sha256")
        if not isinstance(hashes, dict):
            raise ValueError("Unsupported bundle checksum schema")
        for name, expected in hashes.items():
            relative = Path(name)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("Unsafe bundle path")
            if track(directory / relative) != expected:
                raise ValueError(f"Bundle checksum mismatch: {name}")
        return manifest

    base_manifest = verify_bundle(base)
    augmentation_manifest = verify_bundle(augmentation)

    if base_manifest.get("node") != "source_anchored_claim_ledger_v1":
        raise ValueError("Unexpected base ledger")
    if augmentation_manifest.get("node") != (
        "primary_mechanistic_source_collection_v1"
    ):
        raise ValueError("Unexpected augmentation bundle")

    base_sources = read_table(base / "sources.tsv")
    added_sources = read_table(augmentation / "sources.tsv")

    for row in added_sources:
        relative = Path(row["snapshot"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Unsafe augmentation snapshot path")
        path = augmentation / relative
        if track(path) != row["snapshot_sha256"]:
            raise ValueError("Source snapshot hash does not match source row")
        row["snapshot"] = str(path)

    for row in base_sources:
        if track(Path(row["snapshot"])) != row["snapshot_sha256"]:
            raise ValueError("Base source snapshot hash mismatch")

    sources, snapshots = combine_sources(base_sources, added_sources)
    base_claims = read_table(base / "claims.tsv")
    added_claims = read_table(augmentation / "mechanism_claims.tsv")

    for claim in added_claims:
        locator = claim["source_locator"]
        if locator.startswith("sources/"):
            claim["source_locator"] = str(augmentation) + "/" + locator
        if claim["review_status"] != (
            "source_anchors_checked_adjudication_pending"
        ):
            raise ValueError("Unexpected mechanism-claim review status")

    claims = combine_claims(base_claims, added_claims)
    source_ids = {row["source_id"] for row in sources}
    claim_ids = {row["claim_id"] for row in claims}
    if any(row["source_id"] not in source_ids for row in claims):
        raise ValueError("Claim without a registered source")

    links = read_table(base / "claim_hypothesis_links.tsv")
    hypotheses = {}
    identity_fields = [
        "gene", "variant", "selection_witness",
        "q_j_raw", "case_limitations",
    ]
    for row in links:
        identity = {field: row[field] for field in identity_fields}
        hypothesis_id = row["hypothesis_id"]
        if (
            hypothesis_id in hypotheses
            and hypotheses[hypothesis_id] != identity
        ):
            raise ValueError("Inconsistent base hypothesis-link metadata")
        hypotheses[hypothesis_id] = identity

    existing_link_keys = {
        (row["claim_id"], row["hypothesis_id"]) for row in links
    }
    if len(existing_link_keys) != len(links):
        raise ValueError("Duplicate base claim/hypothesis links")

    base_claim_ids = {row["claim_id"] for row in base_claims}
    for claim in added_claims:
        if claim["claim_id"] in base_claim_ids:
            continue
        for hypothesis_id, identity in hypotheses.items():
            links.append({
                "claim_id": claim["claim_id"],
                "hypothesis_id": hypothesis_id,
                **identity,
                "relationship": mechanism_relation(claim["compound"]),
                "applicability_review_status": "pending_review",
            })

    queue = read_table(base / "review_queue.tsv")
    queued = {row["claim_id"] for row in queue}
    for claim in added_claims:
        if claim["claim_id"] not in queued:
            queue.append({
                "claim_id": claim["claim_id"],
                "source_id": claim["source_id"],
                "evidence_class": claim["evidence_class"],
                "required_review": (
                    "experimental system, assay controls, mechanism "
                    "inference and case transferability; "
                    "no treatment nomination from comparator activity"
                ),
                "status": "pending",
            })
            queued.add(claim["claim_id"])

    metadata = load_json(augmentation / "article_metadata.json")
    methods_source_ids = {
        row["source_id"] for row in metadata
        if row["role"] == "rtdetective_primary_methods"
    }

    claim_source_links = [{
        "claim_id": claim["claim_id"],
        "source_id": claim["source_id"],
        "relationship": "primary_source_for_reported_claim",
        "review_status": "scientific_interpretation_pending",
    } for claim in claims]

    for claim in claims:
        if claim["evidence_class"] == "external_precomputed_model_prediction":
            for source_id in sorted(methods_source_ids):
                claim_source_links.append({
                    "claim_id": claim["claim_id"],
                    "source_id": source_id,
                    "relationship": (
                        "methodological_background_not_export_scale_verification"
                    ),
                    "review_status": "pending_methods_review",
                })

    source_review = read_table(augmentation / "review_queue.tsv")

    references = ris_records(
        (base / "used_references.ris").read_text(encoding="utf-8")
    )
    added_references = ris_records(
        (augmentation / "references.ris").read_text(encoding="utf-8")
    )
    for source_id, entry in added_references.items():
        if source_id not in references:
            references[source_id] = entry
    if set(references) != source_ids:
        raise ValueError("Reference IDs and source registry disagree")

    if {row["claim_id"] for row in queue} != claim_ids:
        raise ValueError("Claim review queue is incomplete")
    if len(queue) != len(claim_ids):
        raise ValueError("Duplicate claim review-queue records")
    if any(row["claim_id"] not in claim_ids for row in links):
        raise ValueError("Hypothesis link without a claim")

    link_keys = [
        (row["claim_id"], row["hypothesis_id"]) for row in links
    ]
    if len(link_keys) != len(set(link_keys)):
        raise ValueError("Duplicate merged hypothesis links")

    for row in claim_source_links:
        if (
            row["claim_id"] not in claim_ids
            or row["source_id"] not in source_ids
        ):
            raise ValueError("Invalid claim/source link")

    track(Path(__file__))
    for filename, expected in tracked.items():
        if sha(Path(filename).read_bytes()) != expected:
            raise ValueError("Input changed during augmentation")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    tables = [
        ("sources.tsv", sources, SOURCE_FIELDS),
        ("source_snapshots.tsv", snapshots, SNAPSHOT_FIELDS),
        ("claims.tsv", claims, CLAIM_FIELDS),
        ("claim_hypothesis_links.tsv", links, LINK_FIELDS),
        ("claim_source_links.tsv", claim_source_links, CLAIM_SOURCE_FIELDS),
        ("review_queue.tsv", queue, REVIEW_FIELDS),
        ("source_review_queue.tsv", source_review, SOURCE_REVIEW_FIELDS),
    ]

    with tempfile.TemporaryDirectory(prefix=".merge_", dir=out) as temporary:
        stage = Path(temporary)
        for filename, rows, fields in tables:
            write_table(stage / filename, rows, fields)

        (stage / "used_references.ris").write_text(
            "\n\n".join(
                references[source_id] for source_id in sorted(references)
            ) + "\n",
            encoding="utf-8",
        )
        manifest = {
            "node": "claim_ledger_mechanistic_augmentation_v1",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "input_and_implementation_sha256": tracked,
            "sources": len(sources),
            "source_snapshots": len(snapshots),
            "claims": len(claims),
            "hypothesis_links": len(links),
            "claim_source_links": len(claim_source_links),
            "claim_review_rows": len(queue),
            "source_review_rows": len(source_review),
            "base_ledger_modified": False,
            "automatic_scientific_adjudication": False,
            "comparator_activity_nominates_treatment": False,
            "export_score_units_verified": False,
            "clinical_suitability_assessed": False,
            "weights_modified": False,
            "network_requests": False,
            "output_sha256": {
                path.name: sha(path.read_bytes())
                for path in stage.iterdir() if path.is_file()
            },
        }
        (stage / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )
        for path in list(stage.iterdir()):
            if path.name != "manifest.json":
                path.replace(out / path.name)
        (stage / "manifest.json").replace(out / "manifest.json")

    print(json.dumps({
        "sources": len(sources),
        "claims": len(claims),
        "hypothesis_links": len(links),
        "claim_source_links": len(claim_source_links),
        "claim_review_rows": len(queue),
        "source_review_rows": len(source_review),
    }, sort_keys=True))

if __name__ == "__main__":
    main()
