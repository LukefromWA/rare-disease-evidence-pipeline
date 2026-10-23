#!/usr/bin/env python3
import argparse
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import re
import tempfile

HELPER_PATH = Path(__file__).with_name("d3_readthrough_route.py")
spec = importlib.util.spec_from_file_location(
    "_readthrough_base", HELPER_PATH
)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)

SUPPORTED_VIEWS = {
    "stan", "a10_informed", "a10_agnostic", "dependence_ensemble"
}
WITNESS_FIELDS = [
    "hypothesis_id", "gene", "variant", "selecting_model_view",
    "source_rank", "source_score", "source_score_type",
    "source_field", "source_sha256",
]
COVERAGE_FIELDS = [
    "model_view", "hypotheses", "scored", "not_scored",
    "rank1_hypotheses", "coverage_status",
]

def sha(data):
    return hashlib.sha256(data).hexdigest()

def select_routes(scenarios, views, selected_views, mode="auto"):
    if mode not in {"auto", "off"}:
        raise ValueError("Unsupported routing mode")
    if (
        not isinstance(selected_views, list)
        or not selected_views
        or not all(isinstance(value, str) for value in selected_views)
        or len(selected_views) != len(set(selected_views))
        or not set(selected_views).issubset(SUPPORTED_VIEWS)
    ):
        raise ValueError("Invalid or duplicate model-view selection")

    routed = base.route(scenarios, mode="auto", scope="all")
    hypotheses = {row["hypothesis_id"]: row for row in routed}
    indexed = {}

    required = {
        "hypothesis_id", "gene", "variant", "model_view",
        "view_status", "score_type", "score", "rank",
        "source_field", "source_sha256",
    }

    for row in views:
        if not required.issubset(row):
            raise ValueError("Missing model-view columns")

        if row["model_view"] not in selected_views:
            continue

        key = (row["hypothesis_id"], row["model_view"])
        if key in indexed:
            raise ValueError("Duplicate hypothesis/model-view row")

        hypothesis = hypotheses.get(row["hypothesis_id"])
        if hypothesis is None:
            raise ValueError("Model view refers to an unknown hypothesis")
        if (
            row["gene"] != hypothesis["gene"]
            or row["variant"] != hypothesis["variant"]
        ):
            raise ValueError("Model-view identity mismatch")

        if row["view_status"] == "scored":
            try:
                rank = Decimal(row["rank"])
            except (InvalidOperation, TypeError) as error:
                raise ValueError("Invalid scored-view rank") from error
            if (
                not rank.is_finite()
                or rank < 1
                or rank != rank.to_integral_value()
            ):
                raise ValueError("Scored-view rank must be a positive integer")

            try:
                score = Decimal(row["score"])
            except (InvalidOperation, TypeError) as error:
                raise ValueError("Invalid scored-view score") from error
            if not score.is_finite():
                raise ValueError("Nonfinite scored-view score")

            if (
                not row["score_type"].strip()
                or not row["source_field"].strip()
                or not re.fullmatch(
                    r"[0-9a-fA-F]{64}", row["source_sha256"]
                )
            ):
                raise ValueError("Incomplete scored-view provenance")

        indexed[key] = row

    witnesses = []
    coverage = []

    for view in selected_views:
        missing = [
            hypothesis_id for hypothesis_id in hypotheses
            if (hypothesis_id, view) not in indexed
        ]
        if missing:
            raise ValueError(
                f"Missing hypothesis/model-view records for {view}: "
                f"{len(missing)}"
            )

        records = [
            indexed[(hypothesis_id, view)]
            for hypothesis_id in hypotheses
        ]
        scored = [
            row for row in records if row["view_status"] == "scored"
        ]
        primary = [
            row for row in scored if Decimal(row["rank"]) == 1
        ]

        coverage.append({
            "model_view": view,
            "hypotheses": len(records),
            "scored": len(scored),
            "not_scored": len(records) - len(scored),
            "rank1_hypotheses": len(primary),
            "coverage_status": (
                "no_hypotheses" if not records else
                "no_scored_hypotheses" if not scored else
                "no_rank1_in_handoff" if not primary else
                "rank1_present_with_unscored_hypotheses"
                if len(scored) != len(records) else
                "rank1_present_all_hypotheses_scored"
            ),
        })

        for row in primary:
            witnesses.append({
                "hypothesis_id": row["hypothesis_id"],
                "gene": row["gene"],
                "variant": row["variant"],
                "selecting_model_view": view,
                "source_rank": row["rank"],
                "source_score": row["score"],
                "source_score_type": row["score_type"],
                "source_field": row["source_field"],
                "source_sha256": row["source_sha256"],
            })

    for row in routed:
        selecting_views = [
            view for view in selected_views
            if any(
                witness["hypothesis_id"] == row["hypothesis_id"]
                and witness["selecting_model_view"] == view
                for witness in witnesses
            )
        ]
        selected = bool(selecting_views)

        row["selection_policy"] = "primary_per_view"
        row["selecting_views"] = ";".join(selecting_views)
        row["selected_by_policy"] = selected
        row["selected_for_assessment"] = selected and mode == "auto"

        if mode == "off":
            row["route_status"] = "disabled"
            row["reason"] = "module_disabled_by_configuration"
        elif not selected:
            row["route_status"] = "not_selected"
            row["reason"] = "no_rank1_witness_in_configured_views"

    return routed, witnesses, coverage

def read_table(path):
    raw = path.read_bytes()
    data = gzip.decompress(raw) if path.suffix == ".gz" else raw
    reader = csv.DictReader(
        io.StringIO(data.decode("utf-8-sig")), delimiter="\t"
    )
    if not reader.fieldnames:
        raise ValueError("Missing input table header")
    if len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError("Duplicate input table headers")
    rows = list(reader)
    if any(
        None in row or any(value is None for value in row.values())
        for row in rows
    ):
        raise ValueError("Malformed input table rows")
    return raw, reader.fieldnames, rows

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
    parser.add_argument("--scenarios", required=True)
    parser.add_argument("--model-views", required=True)
    parser.add_argument("--settings-json", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    settings = json.loads(args.settings_json)
    if set(settings) != {"mode", "views"}:
        raise ValueError("Settings must contain exactly mode and views")

    scenario_path = Path(args.scenarios)
    view_path = Path(args.model_views)
    scenario_raw, scenario_fields, scenarios = read_table(scenario_path)
    view_raw, view_fields, views = read_table(view_path)

    if not {
        "hypothesis_id", "gene", "variant", "q_j", "consequence"
    }.issubset(scenario_fields):
        raise ValueError("Incomplete scenario schema")
    if not {
        "hypothesis_id", "gene", "variant", "model_view",
        "view_status", "score_type", "score", "rank",
        "source_field", "source_sha256",
    }.issubset(view_fields):
        raise ValueError("Incomplete model-view schema")

    tracked_paths = [
        scenario_path, view_path, Path(__file__), HELPER_PATH
    ]
    before = {
        str(path): sha(path.read_bytes()) for path in tracked_paths
    }

    if before[str(scenario_path)] != sha(scenario_raw):
        raise ValueError("Scenarios changed while being read")
    if before[str(view_path)] != sha(view_raw):
        raise ValueError("Model views changed while being read")

    routed, witnesses, coverage = select_routes(
        scenarios, views, settings["views"], settings["mode"]
    )
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    route_fields = base.FIELDS + [
        "selection_policy", "selecting_views", "selected_by_policy"
    ]
    tables = [
        ("routing.tsv", routed, route_fields),
        ("selection_witnesses.tsv", witnesses, WITNESS_FIELDS),
        ("view_coverage.tsv", coverage, COVERAGE_FIELDS),
    ]

    with tempfile.TemporaryDirectory(
        prefix=".multiview_", dir=out
    ) as temporary:
        stage = Path(temporary)
        for filename, rows, fields in tables:
            write_table(stage / filename, rows, fields)

        after = {
            str(path): sha(path.read_bytes()) for path in tracked_paths
        }
        if before != after:
            raise ValueError("An input or implementation changed during execution")

        manifest = {
            "node": "readthrough_multiview_routing_v1",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "settings": settings,
            "selection_policy": "primary_per_view",
            "selection_definition": (
                "union of reported rank-1 hypotheses in configured views; "
                "ties retained; absent rank-1 candidates not substituted"
            ),
            "analysis_role": "exploratory_multiview_screening",
            "policy_introduced_after_case_model_disagreement_review": True,
            "input_and_implementation_sha256": before,
            "unique_hypotheses": len(routed),
            "selected_hypotheses": sum(
                row["selected_by_policy"] for row in routed
            ),
            "selection_witnesses": len(witnesses),
            "status_counts": dict(Counter(
                row["route_status"] for row in routed
            )),
            "coverage": coverage,
            "weights_or_scores_modified": False,
            "scores_averaged_or_renormalized": False,
            "genotype_phase_inferred": False,
            "context_verified_by_this_node": False,
            "RTDetective_lookup_performed": False,
            "clinical_eligibility_assessed": False,
            "network_requests": False,
            "output_sha256": {
                filename: sha((stage / filename).read_bytes())
                for filename, _, _ in tables
            },
        }
        (stage / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )

        for filename, _, _ in tables:
            (stage / filename).replace(out / filename)
        (stage / "manifest.json").replace(out / "manifest.json")

    print(json.dumps({
        "selected_hypotheses": manifest["selected_hypotheses"],
        "status_counts": manifest["status_counts"],
    }, sort_keys=True))

if __name__ == "__main__":
    main()
