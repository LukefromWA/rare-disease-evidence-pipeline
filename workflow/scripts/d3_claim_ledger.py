#!/usr/bin/env python3
import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import re
import tempfile
import unicodedata

CLAIM_FIELDS = [
    "claim_id", "source_id", "compound", "evidence_class",
    "source_gene", "population_or_model", "outcome", "claim_text",
    "statistics_json", "source_locator", "source_excerpt",
    "limitations", "dependence_group", "extraction_method", "review_status",
]
SOURCE_FIELDS = [
    "source_id", "source_type", "title", "authors", "year",
    "publication", "doi", "url", "snapshot", "snapshot_sha256",
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

def sha(data):
    return hashlib.sha256(data).hexdigest()

def normalize(value):
    value = unicodedata.normalize("NFKC", str(value))
    for character in ["−", "–", "—"]:
        value = value.replace(character, "-")
    return " ".join(value.split())

def stable_id(prefix, value):
    raw = json.dumps(
        value, sort_keys=True, ensure_ascii=False
    ).encode("utf-8")
    return prefix + sha(raw)[:24]

def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        rows = list(reader)
    if any(
        None in row or any(value is None for value in row.values())
        for row in rows
    ):
        raise ValueError(f"Malformed table: {path}")
    return rows

def write_table(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fields,
            delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

def make_claim(seed, metadata, table_rows):
    pmid = str(seed["pmid"])
    if pmid not in metadata:
        raise ValueError(f"Publication metadata missing: {pmid}")

    anchor = seed["anchor"]
    if anchor["kind"] == "abstract":
        abstract = metadata[pmid].get("abstract", "")
        quotes = anchor["quotes"]
        if not quotes or any(
            normalize(quote) not in normalize(abstract)
            for quote in quotes
        ):
            raise ValueError(
                f"Abstract anchor mismatch: {seed['seed_key']}"
            )
        excerpt = " | ".join(quotes)
        locator = f"publication_metadata.json#/{pmid}/abstract"

    elif anchor["kind"] == "table_cells":
        excerpts = []
        for cell in anchor["cells"]:
            matches = [
                row for row in table_rows
                if row["pmid"] == pmid
                and row["table_id"] == anchor["table_id"]
                and int(row["row_number"]) == cell["row"]
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"Table row missing/ambiguous: {seed['seed_key']}"
                )
            cells = json.loads(matches[0]["cells_json"])
            column = cell["column"]
            if column < 0 or column >= len(cells):
                raise ValueError("Table column outside extracted row")
            actual = cells[column]

            # A separate superscript footnote can have been flattened
            # into a numeric table cell, e.g. "0.8567c".
            if normalize(actual) != normalize(cell["expected"]):
                stripped = re.sub(r"[a-z]$", "", normalize(actual))
                if (
                    not re.fullmatch(
                        r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)",
                        normalize(cell["expected"]),
                    )
                    or stripped != normalize(cell["expected"])
                ):
                    raise ValueError(
                        f"Table cell mismatch: {seed['seed_key']}; "
                        f"row={cell['row']} column={column}; "
                        f"actual={actual!r} expected={cell['expected']!r}"
                    )
            excerpts.append({
                "row": cell["row"],
                "column_zero_based": column,
                "raw_cell": actual,
            })
        excerpt = json.dumps(excerpts, ensure_ascii=False)
        locator = (
            f"fulltext_table_rows.tsv:PMID={pmid};"
            f"table={anchor['table_id']}"
        )
    else:
        raise ValueError("Unsupported evidence anchor kind")

    return {
        "claim_id": stable_id("CLM_", seed),
        "source_id": "PMID" + pmid,
        "compound": seed["compound"],
        "evidence_class": seed["evidence_class"],
        "source_gene": seed["source_gene"],
        "population_or_model": seed["population_or_model"],
        "outcome": seed["outcome"],
        "claim_text": seed["claim_text"],
        "statistics_json": json.dumps(
            seed["statistics"], sort_keys=True
        ),
        "source_locator": locator,
        "source_excerpt": excerpt,
        "limitations": seed["limitations"],
        "dependence_group": "PMID" + pmid,
        "extraction_method": (
            "curated_seed_checked_against_saved_source_anchor"
        ),
        "review_status": "source_anchors_checked_adjudication_pending",
    }


# BEGIN ledger_hypothesis_source_v1_1
def hypotheses_from_bridge(rows):
    hypotheses = {}
    fields = [
        "gene", "variant", "q_j_raw", "selecting_views",
        "phase_status", "haplotype_scope", "NMD_status",
    ]
    for row in rows:
        if row.get("upstream_route_status") != (
            "eligible_for_context_verification"
        ):
            continue
        if (
            not row.get("hypothesis_id")
            or any(field not in row for field in fields)
        ):
            raise ValueError("Incomplete eligible bridge hypothesis")
        hypothesis_id = row["hypothesis_id"]
        if hypothesis_id in hypotheses:
            raise ValueError("Duplicate eligible bridge hypothesis")
        hypotheses[hypothesis_id] = {
            field: row[field] for field in fields
        }
    return hypotheses
# END ledger_hypothesis_source_v1_1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bridge-dir", required=True)
    parser.add_argument("--adoption-dir", required=True)
    parser.add_argument("--fulltext-dir", required=True)
    parser.add_argument("--bibliography-dir", required=True)
    parser.add_argument("--inventory-metadata", required=True)
    parser.add_argument("--seeds", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    adoption = Path(args.adoption_dir)
    fulltext = Path(args.fulltext_dir)
    bibliography = Path(args.bibliography_dir)
    tracked = {}

    def track(path):
        path = Path(path)
        digest = sha(path.read_bytes())
        previous = tracked.get(str(path))
        if previous is not None and previous != digest:
            raise ValueError("Input changed while being read")
        tracked[str(path)] = digest
        return digest

    def load_json(path):
        track(path)
        return json.loads(Path(path).read_text(encoding="utf-8"))

    def verify_bundle(directory):
        manifest = load_json(directory / "manifest.json")
        hashes = manifest.get("output_sha256")
        if not isinstance(hashes, dict):
            raise ValueError("Unsupported bundle checksum schema")
        for filename, expected in hashes.items():
            relative = Path(filename)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("Unsafe bundle output path")
            if track(directory / relative) != expected:
                raise ValueError(f"Bundle checksum mismatch: {filename}")
        return manifest

    adoption_manifest = verify_bundle(adoption)
    if adoption_manifest.get("node") != "rtdetective_cached_lookup_adoption_v1":
        raise ValueError("Unexpected adoption node")
    bridge_dir = Path(args.bridge_dir)
    bridge_manifest = verify_bundle(bridge_dir)
    if bridge_manifest.get("node") != "readthrough_cached_context_bridge_v1":
        raise ValueError("Unexpected bridge node")
    bridge_rows = read_table(bridge_dir / "verification.tsv")
    verify_bundle(fulltext)
    verify_bundle(bibliography)

    metadata = load_json(fulltext / "publication_metadata.json")
    table_rows = read_table(fulltext / "fulltext_table_rows.tsv")
    library = {
        row["citation_key"]: row
        for row in read_table(bibliography / "publication_provenance.tsv")
    }
    seed_bundle = load_json(Path(args.seeds))
    if seed_bundle.get("schema_version") != 1:
        raise ValueError("Unsupported curation-seed schema")
    seeds = seed_bundle["claims"]
    seed_keys = [seed["seed_key"] for seed in seeds]
    if len(seed_keys) != len(set(seed_keys)):
        raise ValueError("Duplicate curation seed keys")

    inventory_path = Path(args.inventory_metadata)
    inventory = load_json(inventory_path)
    inventory_pins = [
        item["sha256"]
        for filename, item in adoption_manifest[
            "input_and_implementation_fingerprints"
        ].items()
        if Path(filename).resolve() == inventory_path.resolve()
    ]
    if inventory_pins != [tracked[str(inventory_path)]]:
        raise ValueError("Inventory is not uniquely pinned by adoption")

    dataset_id = (
        f"FIGSHARE{inventory['id']}v{inventory['version']}"
    )

    claims = [
        make_claim(seed, metadata, table_rows) for seed in seeds
    ]
    source_rows = {}
    for seed in seeds:
        source_id = "PMID" + str(seed["pmid"])
        if source_id not in library:
            raise ValueError(f"Bibliography entry missing: {source_id}")
        row = library[source_id]
        source_rows[source_id] = {
            "source_id": source_id,
            "source_type": "publication",
            "title": row["title"],
            "authors": row["authors"],
            "year": row["year"],
            "publication": row["publication"],
            "doi": row["doi"],
            "url": row["url"],
            "snapshot": str(fulltext / "publication_metadata.json"),
            "snapshot_sha256": tracked[
                str(fulltext / "publication_metadata.json")
            ],
        }

    predictions = read_table(adoption / "predictions.tsv")
    if predictions:
        author_names = []
        for author in inventory.get("authors", []):
            if isinstance(author, dict):
                name = author.get("full_name") or author.get("name")
                if name:
                    author_names.append(name)
        source_rows[dataset_id] = {
            "source_id": dataset_id,
            "source_type": "versioned_prediction_dataset",
            "title": inventory.get("title", ""),
            "authors": "; ".join(author_names),
            "year": "",
            "publication": "Figshare",
            "doi": inventory.get("doi", ""),
            "url": (
                inventory.get("url_public_html")
                or f"https://figshare.com/articles/dataset/{inventory['id']}"
            ),
            "snapshot": str(inventory_path),
            "snapshot_sha256": tracked[str(inventory_path)],
        }

    hypotheses = hypotheses_from_bridge(bridge_rows)
    prediction_claim_links = []
    seen_prediction_claims = set()

    for prediction in predictions:
        hypothesis_id = prediction["hypothesis_id"]
        identity_fields = [
            "gene", "variant", "q_j_raw", "selecting_views",
            "phase_status", "haplotype_scope", "NMD_status",
        ]
        identity = {key: prediction[key] for key in identity_fields}
        if hypothesis_id not in hypotheses:
            raise ValueError("Prediction has no eligible bridge hypothesis")
        if hypotheses[hypothesis_id] != identity:
            raise ValueError("Prediction and bridge hypothesis metadata disagree")

        claim_identity = {
            "source": dataset_id,
            "track": prediction["track_name"],
            "lookup_manifest": prediction["lookup_manifest_sha256"],
            "source_line": prediction["source_line_number"],
            "chrom": prediction["chrom"],
            "start": prediction["target_start"],
            "end": prediction["target_end"],
        }
        claim_id = stable_id("CLM_", claim_identity)

        if claim_id not in seen_prediction_claims:
            claims.append({
                "claim_id": claim_id,
                "source_id": dataset_id,
                "compound": prediction["compound"],
                "evidence_class": "external_precomputed_model_prediction",
                "source_gene": prediction["gene"],
                "population_or_model": "author-exported prediction track",
                "outcome": "raw score for exact cached genomic interval",
                "claim_text": (
                    "Verified cached lookup contains the preserved "
                    "author-exported raw score for the queried interval."
                ),
                "statistics_json": json.dumps({
                    "author_score_raw": prediction["author_score_raw"],
                    "score_units": prediction["score_units"],
                }, sort_keys=True),
                "source_locator": json.dumps(
                    claim_identity, sort_keys=True
                ),
                "source_excerpt": prediction["author_score_raw"],
                "limitations": (
                    "Not a direct measurement of functional rescue. "
                    "Score units unasserted; author-context equivalence "
                    "not independently established."
                ),
                "dependence_group": dataset_id,
                "extraction_method": "verified_cached_lookup_adoption",
                "review_status": (
                    "adoption_verified_interpretation_pending"
                ),
            })
            seen_prediction_claims.add(claim_id)
        prediction_claim_links.append((claim_id, hypothesis_id))

    links = []
    def add_link(claim_id, hypothesis_id, relation):
        hypothesis = hypotheses[hypothesis_id]
        links.append({
            "claim_id": claim_id,
            "hypothesis_id": hypothesis_id,
            "gene": hypothesis["gene"],
            "variant": hypothesis["variant"],
            "relationship": relation,
            "selection_witness": hypothesis["selecting_views"],
            "q_j_raw": hypothesis["q_j_raw"],
            "applicability_review_status": "pending_review",
            "case_limitations": json.dumps({
                key: hypothesis[key] for key in [
                    "phase_status", "haplotype_scope", "NMD_status"
                ]
            }, sort_keys=True),
        })

    for claim_id, hypothesis_id in prediction_claim_links:
        add_link(
            claim_id, hypothesis_id,
            "exact_interval_prediction_not_functional_validation",
        )

    for seed, claim in zip(seeds, claims[:len(seeds)]):
        for hypothesis_id in hypotheses:
            relation = (
                "same_gene_source_variant_and_indication_applicability_unreviewed"
                if seed["source_gene"] == hypotheses[hypothesis_id]["gene"]
                else seed["hypothesis_relation"]
            )
            add_link(claim["claim_id"], hypothesis_id, relation)

    claim_ids = [claim["claim_id"] for claim in claims]
    if len(claim_ids) != len(set(claim_ids)):
        raise ValueError("Duplicate claim IDs")
    if any(claim["source_id"] not in source_rows for claim in claims):
        raise ValueError("Claim without registered source")
    link_keys = [
        (row["claim_id"], row["hypothesis_id"]) for row in links
    ]
    if len(link_keys) != len(set(link_keys)):
        raise ValueError("Duplicate claim/hypothesis links")

    queue = [{
        "claim_id": claim["claim_id"],
        "source_id": claim["source_id"],
        "evidence_class": claim["evidence_class"],
        "required_review": (
            "score scale, model conditions, author-context equivalence "
            "and functional transferability"
            if claim["evidence_class"] == "external_precomputed_model_prediction"
            else "methods, analysis hierarchy, uncertainty, dependence, "
                 "limitations and hypothesis transferability"
        ),
        "status": "pending",
    } for claim in claims]

    ris_path = bibliography / "ataluren_evidence.ris"
    ris_text = ris_path.read_text(encoding="utf-8")
    entries = re.split(r"(?m)(?=^TY  - )", ris_text)
    selected_ris = []
    for source_id in sorted(
        source_id for source_id in source_rows
        if source_id.startswith("PMID")
    ):
        matched = [
            entry for entry in entries
            if re.search(
                r"(?m)^ID  - " + re.escape(source_id) + r"\s*$",
                entry,
            )
        ]
        if len(matched) != 1:
            raise ValueError(f"RIS entry missing/ambiguous: {source_id}")
        selected_ris.append(matched[0].strip())

    if dataset_id in source_rows:
        row = source_rows[dataset_id]
        entry = [
            "TY  - DATA",
            "ID  - " + dataset_id,
            "TI  - " + row["title"],
            "UR  - " + row["url"],
            "N1  - Version " + str(inventory["version"])
                + "; prediction dataset, not clinical evidence",
        ]
        if row["doi"]:
            entry.append("DO  - " + row["doi"])
        entry.append("ER  -")
        selected_ris.append("\n".join(entry))

    track(Path(__file__))
    for filename, expected in tracked.items():
        if sha(Path(filename).read_bytes()) != expected:
            raise ValueError(f"Input changed during ledger build: {filename}")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    tables = [
        ("sources.tsv", list(source_rows.values()), SOURCE_FIELDS),
        ("claims.tsv", claims, CLAIM_FIELDS),
        ("claim_hypothesis_links.tsv", links, LINK_FIELDS),
        ("review_queue.tsv", queue, REVIEW_FIELDS),
    ]
    with tempfile.TemporaryDirectory(prefix=".claims_", dir=out) as temporary:
        stage = Path(temporary)
        for filename, rows, fields in tables:
            write_table(stage / filename, rows, fields)
        (stage / "used_references.ris").write_text(
            "\n\n".join(selected_ris) + "\n", encoding="utf-8"
        )
        manifest = {
            "node": "source_anchored_claim_ledger_v1",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "input_and_implementation_sha256": tracked,
            "sources": len(source_rows),
            "claims": len(claims),
            "hypothesis_links": len(links),
            "hypothesis_source": "nonsense_eligible_bridge_records_not_prediction_availability",
            "link_hypotheses": len(hypotheses),
            "review_queue_rows": len(queue),
            "claim_classes": dict(Counter(
                claim["evidence_class"] for claim in claims
            )),
            "abstract_normalization": (
                "Unicode NFKC, whitespace collapse, dash normalization"
            ),
            "table_columns_are_zero_based": True,
            "one_letter_numeric_table_footnotes_preserved_in_excerpt": True,
            "automatic_evidence_adjudication": False,
            "independent_study_count_inferred": False,
            "numeric_evidence_pooled": False,
            "likelihoods_or_pseudo_counts_generated": False,
            "regulatory_or_access_claims_added": False,
            "clinical_suitability_assessed": False,
            "network_requests": False,
            "original_track_verification_repeated": False,
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
        "sources": len(source_rows),
        "claims": len(claims),
        "hypothesis_links": len(links),
        "review_queue": len(queue),
    }, sort_keys=True))

if __name__ == "__main__":
    main()
