#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import tempfile

NONSENSE_TERMS = {"stop_gained", "nonsense", "nonsense_mutation"}
FIELDS = [
    "hypothesis_id", "gene", "variant", "q_j_raw", "consequence",
    "selected_for_assessment", "nonsense_consequence",
    "route_status", "reason",
]

def sha(data):
    return hashlib.sha256(data).hexdigest()

def consequence_terms(value):
    return {
        term.lower()
        for term in re.split(r"[&,;|\s]+", value.strip())
        if term
    }

def route(rows, mode="auto", scope="primary_q_j"):
    if mode not in {"auto", "off"}:
        raise ValueError("mode must be auto or off")
    if scope not in {"primary_q_j", "all"}:
        raise ValueError("scope must be primary_q_j or all")

    hypotheses = {}
    for row in rows:
        required = ["hypothesis_id", "gene", "variant", "q_j", "consequence"]
        if any(key not in row for key in required):
            raise ValueError("Missing required handoff column")

        identity = tuple(str(row[key]).strip() for key in required[:3])
        if not all(identity):
            raise ValueError("Blank hypothesis identity")

        raw_q = str(row["q_j"]).strip()
        try:
            q = Decimal(raw_q)
        except InvalidOperation as error:
            raise ValueError("Invalid q_j") from error
        if not q.is_finite() or not Decimal(0) <= q <= Decimal(1):
            raise ValueError("q_j must be finite and within [0,1]")

        consequence = str(row["consequence"]).strip()
        item = {
            "identity": identity,
            "q": q,
            "q_raw": raw_q,
            "consequence": consequence,
            "terms": consequence_terms(consequence),
        }

        hypothesis_id = identity[0]
        if hypothesis_id in hypotheses:
            previous = hypotheses[hypothesis_id]
            if (
                previous["identity"] != identity
                or previous["q"] != q
                or previous["terms"] != item["terms"]
            ):
                raise ValueError(
                    "Inconsistent metadata across hypothesis scenarios"
                )
            if previous["q_raw"] != raw_q:
                raise ValueError(
                    "Inconsistent q_j textual representation across scenarios"
                )
        else:
            hypotheses[hypothesis_id] = item

    maximum = max(
        (item["q"] for item in hypotheses.values()), default=None
    )
    output = []

    for hypothesis_id in sorted(hypotheses):
        item = hypotheses[hypothesis_id]
        selected = scope == "all" or item["q"] == maximum
        missing_consequence = (
            not item["consequence"]
            or item["consequence"].lower() in {"na", "nan", "none", "."}
        )
        nonsense = bool(item["terms"] & NONSENSE_TERMS)

        if mode == "off":
            status, reason = "disabled", "module_disabled_by_configuration"
        elif not selected:
            status, reason = (
                "not_selected",
                "outside_configured_hypothesis_scope",
            )
        elif missing_consequence:
            status, reason = (
                "blocked_missing_information",
                "selected_hypothesis_has_no_usable_consequence",
            )
        elif not nonsense:
            status, reason = (
                "not_applicable",
                "selected_hypothesis_has_no_recognized_nonsense_term",
            )
        else:
            status, reason = (
                "eligible_for_context_verification",
                "nonsense_annotation_present_context_not_yet_adjudicated",
            )

        output.append(dict(zip(FIELDS, [
            hypothesis_id,
            item["identity"][1],
            item["identity"][2],
            item["q_raw"],
            item["consequence"],
            selected,
            nonsense,
            status,
            reason,
        ])))

    return output

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scenarios", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--mode", choices=["auto", "off"], default="auto")
    parser.add_argument(
        "--scope", choices=["primary_q_j", "all"], default="primary_q_j"
    )
    args = parser.parse_args()

    source = Path(args.scenarios)
    raw = source.read_bytes()
    decoded = gzip.decompress(raw) if source.suffix == ".gz" else raw

    reader = csv.DictReader(
        io.StringIO(decoded.decode("utf-8-sig")), delimiter="\t"
    )
    required = {"hypothesis_id", "gene", "variant", "q_j", "consequence"}
    if not required.issubset(reader.fieldnames or []):
        raise ValueError("Scenarios table lacks required routing columns")

    rows = list(reader)
    output = route(rows, args.mode, args.scope)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".route_", dir=out) as temporary:
        stage = Path(temporary)
        table = stage / "routing.tsv"

        with table.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=FIELDS,
                delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(output)

        settings = {
            "mode": args.mode,
            "scope": args.scope,
            "primary_selection_definition": (
                "maximum original handoff q_j; all ties retained"
            ),
        }
        manifest = {
            "node": "readthrough_annotation_routing_v1",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "input": str(source),
            "input_sha256": sha(raw),
            "implementation_sha256": sha(Path(__file__).read_bytes()),
            "settings": settings,
            "source_scenario_rows": len(rows),
            "unique_hypotheses": len(output),
            "status_counts": dict(Counter(
                row["route_status"] for row in output
            )),
            "hypothesis_weights_modified": False,
            "phase_or_multiallelic_membership_inferred": False,
            "transcript_or_genomic_context_verified_by_this_node": False,
            "RTDetective_lookup_performed": False,
            "network_requests": False,
            "clinical_eligibility_assessed": False,
            "output_sha256": {"routing.tsv": sha(table.read_bytes())},
        }

        if sha(source.read_bytes()) != sha(raw):
            raise ValueError("Scenarios changed during execution")

        (stage / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        table.replace(out / "routing.tsv")
        (stage / "manifest.json").replace(out / "manifest.json")

    print(json.dumps(manifest["status_counts"], sort_keys=True))

if __name__ == "__main__":
    main()
