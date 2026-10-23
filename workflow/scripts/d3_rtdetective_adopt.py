#!/usr/bin/env python3
import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile

CATALOG = {
    "DAP_tga.txt": "2,6-diaminopurine",
    "Clitocine_tga.txt": "clitocine",
    "SRI_tga.txt": "SRI-41315",
}
STATUS_FIELDS = [
    "hypothesis_id", "gene", "variant", "compound", "track_name",
    "adoption_status", "reason",
]
PREDICTION_FIELDS = [
    "hypothesis_id", "gene", "variant", "q_j_raw", "selecting_views",
    "transcript_version", "alternate_codon_RNA",
    "alternate_context_bracketed_RNA",
    "assembly", "chrom", "target_start", "target_end",
    "compound", "track_name", "author_score_raw", "source_line_number",
    "article_id", "article_version",
    "source_track_md5", "source_track_sha256",
    "lookup_manifest_sha256",
    "evidence_type", "score_units",
    "author_context_equivalence",
    "phase_status", "haplotype_scope", "NMD_status",
    "functional_rescue", "clinical_suitability",
]

def fingerprint(path):
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
            size += len(block)
    return {
        "bytes": size,
        "md5": md5.hexdigest(),
        "sha256": sha256.hexdigest(),
    }

def read_table(path):
    with path.open(
        encoding="utf-8-sig", newline=""
    ) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing table header: {path}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Duplicate headers: {path}")
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

def exact_record(bridge, statuses, exact):
    key = (
        bridge["chrom"],
        int(bridge["bed_chromStart"]),
        int(bridge["bed_chromEnd"]),
    )
    matched_statuses = [
        row for row in statuses
        if (
            row["chrom"], int(row["target_start"]),
            int(row["target_end"])
        ) == key
    ]
    if not matched_statuses:
        return None, "no_cached_target_status"
    if len(matched_statuses) != 1:
        return None, "ambiguous_cached_target_status"

    status = matched_statuses[0]
    if (
        status["status"] != "unique_exact_ascending_coordinate_match"
        or int(status["n_exact_records"]) != 1
        or int(status["n_nonforward_records_near_target"]) != 0
    ):
        return None, "cached_lookup_not_unique_and_unambiguous"

    matches = [
        row for row in exact
        if (
            row["chrom"], int(row["target_start"]),
            int(row["target_end"])
        ) == key
    ]
    if len(matches) != 1:
        return None, "exact_record_count_disagrees_with_status"

    record = matches[0]
    if (
        record["coordinate_order"] != "ascending"
        or int(record["raw_start"]) != key[1]
        or int(record["raw_end"]) != key[2]
        or key[2] - key[1] != 3
    ):
        return None, "exact_record_interval_disagreement"

    try:
        finite = math.isfinite(float(record["author_score_raw"]))
        source_line = int(record["source_line_number"])
    except (KeyError, ValueError, OverflowError):
        return None, "invalid_score_or_source_line"
    if not finite or source_line < 1:
        return None, "invalid_score_or_source_line"

    return record, "accepted_unique_exact_record"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bridge-dir", required=True)
    parser.add_argument("--mapping-table", required=True)
    parser.add_argument("--lookup-dir", required=True)
    parser.add_argument("--inventory-dir", required=True)
    parser.add_argument("--track-cache", required=True)
    parser.add_argument("--declared-inputs-json", required=True)
    parser.add_argument("--settings-json", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    settings = json.loads(args.settings_json)
    if set(settings) != {"mode", "tracks", "article_id", "article_version"}:
        raise ValueError("Unexpected adoption settings")
    if settings["mode"] not in {"auto", "off"}:
        raise ValueError("Unsupported adoption mode")
    tracks = settings["tracks"]
    if (
        not isinstance(tracks, list) or not tracks
        or not all(isinstance(name, str) for name in tracks)
        or len(tracks) != len(set(tracks))
        or not set(tracks).issubset(CATALOG)
    ):
        raise ValueError("Unsupported or duplicate track selection")

    declared = {
        str(Path(value)) for value in json.loads(args.declared_inputs_json)
    }
    fingerprints = {}

    def tracked(path):
        path = Path(path)
        name = str(path)
        if name not in declared:
            raise ValueError(f"Undeclared input: {path}")
        if name not in fingerprints:
            fingerprints[name] = fingerprint(path)
        return fingerprints[name]

    def load_json(path):
        tracked(path)
        return json.loads(Path(path).read_text(encoding="utf-8"))

    def verify_outputs(directory, manifest):
        hashes = manifest.get("output_sha256")
        if not isinstance(hashes, dict):
            raise ValueError("Unsupported output-checksum manifest schema")
        for filename, expected in hashes.items():
            relative = Path(filename)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("Unsafe manifest output path")
            if tracked(directory / relative)["sha256"] != expected:
                raise ValueError(f"Output checksum mismatch: {filename}")

    bridge_dir = Path(args.bridge_dir)
    bridge_manifest = load_json(bridge_dir / "manifest.json")
    verify_outputs(bridge_dir, bridge_manifest)

    if bridge_manifest.get("node") != "readthrough_cached_context_bridge_v1":
        raise ValueError("Unexpected bridge node")

    for filename, expected in bridge_manifest[
        "input_and_implementation_sha256"
    ].items():
        path = Path(filename)
        actual = fingerprint(path)
        if actual["sha256"] != expected:
            raise ValueError(
                f"Bridge upstream fingerprint no longer matches: {filename}"
            )
        fingerprints[str(path)] = actual

    bridge_rows = read_table(bridge_dir / "verification.tsv")
    eligible = [
        row for row in bridge_rows
        if row["bridge_status"] == "verified_cached_context_and_mapping"
    ]
    ids = [row["hypothesis_id"] for row in eligible]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate eligible bridge hypotheses")

    predictions, statuses = [], []
    inventory_dir = Path(args.inventory_dir)
    lookup_root = Path(args.lookup_dir)
    track_cache = Path(args.track_cache)
    mapping_path = Path(args.mapping_table)

    metadata = None
    mapping_hash = None

    if eligible and settings["mode"] == "auto":
        mapping_hash = tracked(mapping_path)["sha256"]
        pinned_mapping = [
            expected
            for filename, expected in bridge_manifest[
                "input_and_implementation_sha256"
            ].items()
            if Path(filename).resolve() == mapping_path.resolve()
        ]
        if pinned_mapping != [mapping_hash]:
            raise ValueError("Mapping table is not uniquely pinned by bridge")

        metadata_path = inventory_dir / "article_metadata.json"
        inventory_manifest_path = inventory_dir / "manifest.json"

        if (
            str(metadata_path) in declared
            and str(inventory_manifest_path) in declared
        ):
            inventory_manifest = load_json(inventory_manifest_path)
            if (
                tracked(metadata_path)["sha256"]
                != inventory_manifest["metadata_sha256"]
            ):
                raise ValueError("Pinned inventory metadata checksum mismatch")
            metadata = load_json(metadata_path)
            if (
                metadata.get("id") != settings["article_id"]
                or metadata.get("version") != settings["article_version"]
            ):
                raise ValueError("Inventory article/version mismatch")

    loaded = {}

    def get_lookup(filename):
        if filename in loaded:
            return loaded[filename]

        directory = lookup_root / Path(filename).stem
        manifest_path = directory / "manifest.json"
        track_path = track_cache / filename

        if str(manifest_path) not in declared:
            result = (None, "missing_cached_lookup")
        elif metadata is None:
            result = (None, "missing_pinned_inventory")
        elif str(track_path) not in declared:
            result = (None, "missing_cached_track")
        else:
            entries = [
                item for item in metadata["files"]
                if item["name"] == filename
            ]
            if len(entries) != 1:
                raise ValueError("Track not uniquely identified in inventory")
            entry = entries[0]
            published_md5 = (
                entry.get("computed_md5") or entry.get("supplied_md5") or ""
            ).lower()
            if len(published_md5) != 32:
                raise ValueError("Missing published track MD5")

            print(f"Hash-verifying cached track: {filename}", flush=True)
            actual = tracked(track_path)
            if (
                actual["bytes"] != int(entry["size"])
                or actual["md5"] != published_md5
            ):
                raise ValueError("Cached track fails published size/MD5")

            manifest = load_json(manifest_path)
            verify_outputs(directory, manifest)

            if (
                manifest.get("node")
                != "rtdetective_direction_preserving_coordinate_lookup"
                or manifest.get("source_article_id") != settings["article_id"]
                or manifest.get("source_article_version")
                != settings["article_version"]
                or manifest.get("mapping_table_sha256") != mapping_hash
                or manifest.get("published_md5_verified") is not True
                or manifest.get("source_track_md5") != actual["md5"]
                or manifest.get("source_track_sha256") != actual["sha256"]
                or manifest.get(
                    "descending_intervals_accepted_as_forward_matches"
                ) is not False
                or manifest.get("duplicate_matches_aggregated") is not False
            ):
                raise ValueError(
                    f"Lookup provenance/integrity policy mismatch: {filename}"
                )

            recorded_track = Path(manifest["source_track"])
            if recorded_track.resolve() != track_path.resolve():
                raise ValueError("Lookup source-track path identity mismatch")

            for required in ["lookup_status.tsv", "exact_matches.tsv"]:
                if required not in manifest["output_sha256"]:
                    raise ValueError("Required lookup table is not checksum-pinned")

            result = ({
                "manifest": manifest,
                "manifest_sha256": tracked(manifest_path)["sha256"],
                "track": actual,
                "statuses": read_table(directory / "lookup_status.tsv"),
                "exact": read_table(directory / "exact_matches.tsv"),
            }, "")
        loaded[filename] = result
        return result

    for bridge in eligible:
        for filename in tracks:
            status = {
                "hypothesis_id": bridge["hypothesis_id"],
                "gene": bridge["gene"],
                "variant": bridge["variant"],
                "compound": CATALOG[filename],
                "track_name": filename,
                "adoption_status": "",
                "reason": "",
            }

            if settings["mode"] == "off":
                status.update(
                    adoption_status="disabled",
                    reason="adoption_disabled_by_configuration",
                )
            elif bridge["alternate_codon_RNA"] != "UGA":
                status.update(
                    adoption_status="unsupported_resource_scope",
                    reason="configured_tracks_are_UGA_specific",
                )
            else:
                lookup, problem = get_lookup(filename)
                if lookup is None:
                    status.update(
                        adoption_status=problem,
                        reason="no_prediction_assigned",
                    )
                else:
                    record, reason = exact_record(
                        bridge, lookup["statuses"], lookup["exact"]
                    )
                    if record is None:
                        status.update(
                            adoption_status="lookup_requires_review",
                            reason=reason,
                        )
                    else:
                        prediction = {
                            "hypothesis_id": bridge["hypothesis_id"],
                            "gene": bridge["gene"],
                            "variant": bridge["variant"],
                            "q_j_raw": bridge["q_j_raw"],
                            "selecting_views": bridge["selecting_views"],
                            "transcript_version": bridge["transcript_version"],
                            "alternate_codon_RNA": bridge["alternate_codon_RNA"],
                            "alternate_context_bracketed_RNA":
                                bridge["alternate_context_bracketed_RNA"],
                            "assembly": bridge["assembly"],
                            "chrom": bridge["chrom"],
                            "target_start": bridge["bed_chromStart"],
                            "target_end": bridge["bed_chromEnd"],
                            "compound": CATALOG[filename],
                            "track_name": filename,
                            "author_score_raw": record["author_score_raw"],
                            "source_line_number": record["source_line_number"],
                            "article_id": settings["article_id"],
                            "article_version": settings["article_version"],
                            "source_track_md5": lookup["track"]["md5"],
                            "source_track_sha256": lookup["track"]["sha256"],
                            "lookup_manifest_sha256":
                                lookup["manifest_sha256"],
                            "evidence_type":
                                "external_precomputed_model_prediction",
                            "score_units": "not_asserted",
                            "author_context_equivalence":
                                "not_independently_established",
                            "phase_status": bridge["phase_status"],
                            "haplotype_scope": bridge["haplotype_scope"],
                            "NMD_status": bridge["NMD_status"],
                            "functional_rescue": "not_established",
                            "clinical_suitability": "not_assessed",
                        }
                        predictions.append(prediction)
                        status.update(
                            adoption_status="adopted_verified_cached_lookup",
                            reason=reason,
                        )
            statuses.append(status)

    implementation = Path(__file__)
    fingerprints[str(implementation)] = fingerprint(implementation)

    print("Rechecking tracked input fingerprints before publication...", flush=True)
    for filename, expected in fingerprints.items():
        if fingerprint(Path(filename)) != expected:
            raise ValueError(f"Input changed during adoption: {filename}")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    tables = [
        ("predictions.tsv", predictions, PREDICTION_FIELDS),
        ("adoption_status.tsv", statuses, STATUS_FIELDS),
    ]
    with tempfile.TemporaryDirectory(prefix=".adopt_", dir=out) as temporary:
        stage = Path(temporary)
        for filename, rows, fields in tables:
            write_table(stage / filename, rows, fields)
        manifest = {
            "node": "rtdetective_cached_lookup_adoption_v1",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "settings": settings,
            "eligible_bridge_hypotheses": len(eligible),
            "prediction_rows": len(predictions),
            "status_counts": dict(Counter(
                row["adoption_status"] for row in statuses
            )),
            "input_and_implementation_fingerprints": fingerprints,
            "independent_coordinate_rescan_performed": False,
            "score_units_inferred": False,
            "scores_converted_or_averaged": False,
            "hypothesis_weights_modified": False,
            "author_context_equivalence_independently_established": False,
            "network_requests": False,
            "treatment_ranking_performed": False,
            "output_sha256": {
                filename: fingerprint(stage / filename)["sha256"]
                for filename, _, _ in tables
            },
        }
        (stage / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )
        for filename, _, _ in tables:
            (stage / filename).replace(out / filename)
        (stage / "manifest.json").replace(out / "manifest.json")

    print(json.dumps({
        "prediction_rows": len(predictions),
        "status_counts": manifest["status_counts"],
    }, sort_keys=True))

if __name__ == "__main__":
    main()
