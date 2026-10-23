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

FIELDS = [
    "hypothesis_id", "gene", "variant", "q_j_raw", "selecting_views",
    "upstream_route_status", "bridge_status", "reason",
    "transcript_version", "alternate_codon_RNA",
    "alternate_context_bracketed_RNA",
    "assembly", "chrom", "strand", "bed_chromStart", "bed_chromEnd",
    "phase_status", "haplotype_scope", "NMD_status",
]
REQUEST_FIELDS = [
    "hypothesis_id", "gene", "variant", "selecting_views",
    "bridge_status", "reason",
]

def sha(data):
    return hashlib.sha256(data).hexdigest()

def truth(value):
    return str(value).strip().lower() == "true"

def read_table(path):
    raw = path.read_bytes()
    reader = csv.DictReader(
        raw.decode("utf-8-sig").splitlines(), delimiter="\t"
    )
    if not reader.fieldnames:
        raise ValueError(f"Missing table header: {path}")
    if len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError(f"Duplicate table headers: {path}")
    rows = list(reader)
    if any(
        None in row or any(value is None for value in row.values())
        for row in rows
    ):
        raise ValueError(f"Malformed rows: {path}")
    return rows

def check_record(request, contexts, mappings, assembly="GRCh38"):
    result = {field: "" for field in FIELDS}
    result.update({
        "hypothesis_id": request["hypothesis_id"],
        "gene": request["gene"],
        "variant": request["variant"],
        "q_j_raw": request["q_j_raw"],
        "selecting_views": request["selecting_views"],
        "upstream_route_status": request["route_status"],
    })

    def finish(status, reason):
        result["bridge_status"] = status
        result["reason"] = reason
        return result

    if not truth(request["selected_for_assessment"]):
        return finish(
            "not_requested_upstream",
            "hypothesis_not_enabled_for_assessment",
        )
    if request["route_status"] != "eligible_for_context_verification":
        return finish(
            "not_applicable_upstream",
            "upstream_route_does_not_request_this_module",
        )

    key = (request["gene"], request["variant"])
    context_matches = [
        row for row in contexts
        if (row.get("gene"), row.get("candidate_id")) == key
    ]
    if not context_matches:
        return finish("missing_context", "no_exact_gene_and_allele_match")
    if len(context_matches) != 1:
        return finish("ambiguous_context", "multiple_exact_context_records")

    context = context_matches[0]
    for field in [
        "transcript_version", "alternate_codon_RNA",
        "alternate_context_bracketed_RNA",
        "phase_status", "haplotype_scope", "NMD_status",
    ]:
        result[field] = context.get(field, "")

    required_flags = [
        "sequence_verified", "alternate_is_stop",
        "context_boundary_truncated",
    ]
    if any(field not in context for field in required_flags):
        return finish(
            "blocked_missing_information", "required_context_flags_missing"
        )
    if not truth(context["sequence_verified"]):
        return finish("verification_failed", "upstream_sequence_not_verified")
    if not truth(context["alternate_is_stop"]):
        return finish("verification_failed", "alternate_codon_not_verified_stop")
    if str(context["context_boundary_truncated"]).lower() != "false":
        return finish(
            "verification_failed", "context_truncated_or_boundary_flag_unknown"
        )

    stop = context.get("alternate_codon_RNA", "")
    bracketed = context.get("alternate_context_bracketed_RNA", "")
    if stop not in {"UGA", "UAG", "UAA"}:
        return finish("verification_failed", "unrecognized_alternate_stop")
    if not re.fullmatch(r"[ACGU]+\[" + stop + r"\][ACGU]+", bracketed):
        return finish("verification_failed", "invalid_or_incomplete_RNA_context")
    if context.get("reference_assembly") != assembly:
        return finish("unsupported_resource_scope", "context_assembly_mismatch")

    mapping_matches = [
        row for row in mappings
        if (row.get("gene"), row.get("candidate_id")) == key
    ]
    if not mapping_matches:
        return finish("missing_mapping", "no_exact_gene_and_allele_mapping")
    if len(mapping_matches) != 1:
        return finish("ambiguous_mapping", "multiple_exact_mapping_records")

    mapping = mapping_matches[0]
    for field in [
        "assembly", "chrom", "strand", "bed_chromStart", "bed_chromEnd"
    ]:
        result[field] = mapping.get(field, "")

    for field in [
        "transcript_version", "HGVSc", "HGVSp",
        "alternate_codon_RNA", "alternate_context_bracketed_RNA",
    ]:
        if (
            not context.get(field)
            or context.get(field) != mapping.get(field)
        ):
            return finish(
                "identity_or_transcript_mismatch",
                "context_mapping_disagreement:" + field,
            )

    if mapping.get("assembly") != assembly:
        return finish("unsupported_resource_scope", "mapping_assembly_mismatch")
    if mapping.get("strand") != "1":
        return finish(
            "unsupported_resource_scope",
            "this_cached_lookup_adapter_supports_positive_strand_only",
        )

    for field in [
        "genomic_coordinate_mapping_verified",
        "genomic_context_matches_frozen_CDS",
        "codon_contiguous",
    ]:
        if not truth(mapping.get(field, "")):
            return finish(
                "verification_failed", "mapping_check_not_true:" + field
            )
    if mapping.get("track_lookup_status") != "contiguous_stop_codon_lookup_key":
        return finish(
            "unsupported_resource_scope", "mapping_not_lookup_compatible"
        )

    try:
        start = int(mapping["bed_chromStart"])
        end = int(mapping["bed_chromEnd"])
    except (KeyError, ValueError, TypeError):
        return finish("verification_failed", "invalid_interval")

    if start < 0 or end - start != 3:
        return finish("verification_failed", "not_a_forward_three_base_interval")

    variant = request["variant"].split(":")
    if len(variant) != 4:
        return finish(
            "unsupported_resource_scope",
            "variant_key_not_CHROM_POS_REF_ALT",
        )
    chrom, position, ref, alt = variant
    if len(ref) != 1 or len(alt) != 1:
        return finish(
            "unsupported_resource_scope", "this_adapter_requires_an_SNV"
        )
    try:
        position = int(position)
    except ValueError:
        return finish("verification_failed", "invalid_variant_position")

    if (
        chrom.removeprefix("chr")
        != mapping["chrom"].removeprefix("chr")
        or not start <= position - 1 < end
    ):
        return finish(
            "identity_or_transcript_mismatch",
            "variant_position_not_within_mapped_codon",
        )

    return finish(
        "verified_cached_context_and_mapping",
        "exact_allele_binding_and_upstream_verification_flags_pass",
    )

def verify_outputs(directory, manifest):
    output_hashes = manifest.get("output_sha256")
    if not isinstance(output_hashes, dict):
        raise ValueError(
            f"Unsupported manifest schema in {directory}: "
            "expected output_sha256 mapping"
        )
    for name, expected in output_hashes.items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Unsafe output path in manifest")
        path = directory / relative
        if (
            not isinstance(expected, str)
            or not re.fullmatch(r"[0-9a-fA-F]{64}", expected)
            or not path.is_file()
            or sha(path.read_bytes()) != expected.lower()
        ):
            raise ValueError(f"Output integrity check failed: {path}")

def verify_context_provenance(directory, manifest):
    for key in ["input_sha256", "reference_manifests"]:
        entries = manifest.get(key)
        if not isinstance(entries, dict) or not entries:
            raise ValueError(f"Missing context provenance map: {key}")
        for filename, expected in entries.items():
            path = Path(filename)
            if (
                not isinstance(expected, str)
                or not re.fullmatch(r"[0-9a-fA-F]{64}", expected)
                or not path.is_file()
                or sha(path.read_bytes()) != expected.lower()
            ):
                raise ValueError(
                    f"Context upstream provenance check failed: {filename}"
                )

def write_table(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fields,
            delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--route-dir", required=True)
    parser.add_argument("--context-dir", required=True)
    parser.add_argument("--cache-inputs-json", required=True)
    parser.add_argument("--assembly", default="GRCh38")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    route_dir = Path(args.route_dir)
    context_dir = Path(args.context_dir)
    mapping_dir = context_dir / "genomic_mapping"
    declared = {str(Path(value)) for value in json.loads(args.cache_inputs_json)}

    expected_cache = {
        str(context_dir / "codon_context.tsv"),
        str(context_dir / "manifest.json"),
        str(mapping_dir / "codon_genomic_mapping.tsv"),
        str(mapping_dir / "manifest.json"),
    }
    if not declared.issubset(expected_cache):
        raise ValueError("Unexpected declared cache input")

    tracked = [
        route_dir / "routing.tsv",
        route_dir / "selection_witnesses.tsv",
        route_dir / "view_coverage.tsv",
        route_dir / "manifest.json",
        Path(__file__),
    ] + [Path(value) for value in sorted(declared)]
    before = {str(path): sha(path.read_bytes()) for path in tracked}

    route_manifest = json.loads(
        (route_dir / "manifest.json").read_text()
    )
    verify_outputs(route_dir, route_manifest)
    for filename, expected in route_manifest[
        "input_and_implementation_sha256"
    ].items():
        path = Path(filename)
        if not path.is_file() or sha(path.read_bytes()) != expected:
            raise ValueError(
                "Routing upstream input or implementation no longer matches"
            )

    routing = read_table(route_dir / "routing.tsv")
    requests = [
        row for row in routing
        if truth(row["selected_for_assessment"])
        and row["route_status"] == "eligible_for_context_verification"
    ]

    contexts, mappings = [], []
    context_pair = {
        str(context_dir / "codon_context.tsv"),
        str(context_dir / "manifest.json"),
    }
    mapping_pair = {
        str(mapping_dir / "codon_genomic_mapping.tsv"),
        str(mapping_dir / "manifest.json"),
    }

    if requests and context_pair.issubset(declared):
        context_manifest = json.loads(
            (context_dir / "manifest.json").read_text()
        )
        verify_outputs(context_dir, context_manifest)
        if "codon_context.tsv" not in context_manifest["output_sha256"]:
            raise ValueError("Context table lacks an explicit output checksum")
        verify_context_provenance(context_dir, context_manifest)
        contexts = read_table(context_dir / "codon_context.tsv")

    if requests and mapping_pair.issubset(declared):
        mapping_manifest = json.loads(
            (mapping_dir / "manifest.json").read_text()
        )
        verify_outputs(mapping_dir, mapping_manifest)
        if "codon_genomic_mapping.tsv" not in mapping_manifest["output_sha256"]:
            raise ValueError("Mapping table lacks an explicit output checksum")
        mappings = read_table(mapping_dir / "codon_genomic_mapping.tsv")

    selected = [
        row for row in routing if truth(row["selected_by_policy"])
    ]
    results = [
        check_record(row, contexts, mappings, args.assembly)
        for row in selected
    ]
    request_ids = {row["hypothesis_id"] for row in requests}
    request_rows = [
        {field: row[field] for field in REQUEST_FIELDS}
        for row in results if row["hypothesis_id"] in request_ids
    ]

    after = {str(path): sha(path.read_bytes()) for path in tracked}
    if before != after:
        raise ValueError("A declared input changed during execution")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    tables = [
        ("verification.tsv", results, FIELDS),
        ("context_requests.tsv", request_rows, REQUEST_FIELDS),
    ]
    with tempfile.TemporaryDirectory(prefix=".bridge_", dir=out) as temporary:
        stage = Path(temporary)
        for filename, rows, fields in tables:
            write_table(stage / filename, rows, fields)

        manifest = {
            "node": "readthrough_cached_context_bridge_v1",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "input_and_implementation_sha256": before,
            "requested_assembly": args.assembly,
            "selected_hypotheses": len(selected),
            "context_requests": len(request_rows),
            "status_counts": dict(Counter(
                row["bridge_status"] for row in results
            )),
            "verification_scope": (
                "cached allele/context/mapping binding, output checksums, "
                "and upstream context input/reference-manifest fingerprints"
            ),
            "sequence_recomputed_independently": False,
            "upstream_verification_flags_used": True,
            "phase_or_NMD_status_inferred": False,
            "weights_modified": False,
            "network_requests": False,
            "RTDetective_lookup_performed": False,
            "clinical_eligibility_assessed": False,
            "output_sha256": {
                filename: sha((stage / filename).read_bytes())
                for filename, _, _ in tables
            },
        }
        (stage / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )
        for filename, _, _ in tables:
            (stage / filename).replace(out / filename)
        (stage / "manifest.json").replace(out / "manifest.json")

    print(json.dumps(manifest["status_counts"], sort_keys=True))

if __name__ == "__main__":
    main()
