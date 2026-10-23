#!/usr/bin/env python3
"""Verify selected transcript CDS variants and export focal-variant contexts."""
import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import tempfile
import time
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

BASES = "TCAG"
AMINO_ACIDS = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
CODE = dict(zip(
    ("".join(x) for x in itertools.product(BASES, repeat=3)),
    AMINO_ACIDS,
    strict=True,
))
AA3 = dict(zip(
    (
        "Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys "
        "Met Phe Pro Ser Thr Trp Tyr Val Ter"
    ).split(),
    "ARNDCQEGHILKMFPSTWYV*",
    strict=True,
))

def sha(data):
    return hashlib.sha256(data).hexdigest()

def read_tsv(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))

def parse_hgvs(hgvsc, hgvsp):
    c = re.fullmatch(
        r"(ENST[0-9]+)\.([0-9]+):c\.([0-9]+)([ACGT])>([ACGT])",
        hgvsc,
    )
    p = re.fullmatch(
        r"(?:[^:]+:)?p\.\(?([A-Z][a-z]{2})([0-9]+)"
        r"([A-Z][a-z]{2}|\*)\)?",
        hgvsp,
    )
    if not c or not p:
        raise ValueError(
            f"Unsupported HGVS for this SNV context module: {hgvsc}; {hgvsp}"
        )
    tx, version, position, ref, alt = c.groups()
    ref_aa, aa_position, alt_aa = p.groups()
    if ref_aa not in AA3 or (alt_aa != "*" and alt_aa not in AA3):
        raise ValueError("Unsupported protein amino-acid notation")

    return {
        "transcript": tx,
        "version": int(version),
        "position": int(position),
        "ref": ref,
        "alt": alt,
        "aa_position": int(aa_position),
        "ref_aa": AA3[ref_aa],
        "alt_aa": "*" if alt_aa == "*" else AA3[alt_aa],
    }

def validate_lookup(meta, transcript, version):
    if (
        meta.get("id") != transcript
        or int(meta.get("version", -1)) != version
        or meta.get("object_type") != "Transcript"
        or meta.get("species") != "homo_sapiens"
    ):
        raise ValueError(
            "Exact human transcript/version mismatch: "
            f"requested {transcript}.{version}; "
            f"returned {meta.get('id')}.{meta.get('version')}. "
            "Use a matching reference/archive; no substitution performed."
        )

def annotate(cds, hgvsc, hgvsp, flank):
    v = parse_hgvs(hgvsc, hgvsp)
    position = v["position"]
    start = 3 * (v["aa_position"] - 1)

    if (position - 1) // 3 + 1 != v["aa_position"]:
        raise ValueError("Coding position and protein position disagree")
    if not (1 <= position <= len(cds)) or start + 3 > len(cds):
        raise ValueError("Variant outside retrieved CDS")
    if cds[position - 1] != v["ref"]:
        raise ValueError("HGVS reference nucleotide does not match CDS")

    ref_codon = cds[start:start + 3]
    alt_codon = list(ref_codon)
    alt_codon[(position - 1) - start] = v["alt"]
    alt_codon = "".join(alt_codon)

    if CODE.get(ref_codon) != v["ref_aa"]:
        raise ValueError("Reference codon does not match reported amino acid")
    if CODE.get(alt_codon) != v["alt_aa"]:
        raise ValueError("Altered codon does not match reported protein change")

    left = cds[max(0, start - flank):start]
    right = cds[start + 3:start + 3 + flank]
    if set(left + ref_codon + right) - set("ACGT"):
        raise ValueError("Ambiguous bases in requested coding context")

    rna = lambda text: text.replace("T", "U")
    stop = CODE[alt_codon] == "*"

    return {
        "transcript_version": f'{v["transcript"]}.{v["version"]}',
        "HGVSc": hgvsc,
        "HGVSp": hgvsp,
        "coding_position": position,
        "protein_position": v["aa_position"],
        "position_within_codon": (position - 1) % 3 + 1,
        "reference_codon_DNA": ref_codon,
        "alternate_codon_DNA": alt_codon,
        "reference_codon_RNA": rna(ref_codon),
        "alternate_codon_RNA": rna(alt_codon),
        "alternate_is_stop": stop,
        "upstream_bases": len(left),
        "downstream_bases": len(right),
        "reference_context_DNA": left + ref_codon + right,
        "alternate_context_DNA": left + alt_codon + right,
        "reference_context_RNA": rna(left + ref_codon + right),
        "alternate_context_RNA": rna(left + alt_codon + right),
        "alternate_context_bracketed_RNA":
            rna(left) + "[" + rna(alt_codon) + "]" + rna(right),
        "stop_plus4_base_RNA": rna(right[:1]) if stop else "",
        "stop_plus4_motif_RNA":
            rna(alt_codon + right[:1]) if stop and right else "",
        "context_boundary_truncated":
            len(left) < flank or len(right) < flank,
        "sequence_verified": True,
        "NMD_status": "not_assessed",
        "genomic_coordinate_mapping_verified": False,
        "haplotype_scope": "focal_variant_only_other_variants_not_applied",
    }

def fetch_json(server, endpoint, params=None):
    url = server.rstrip("/") + endpoint
    if params:
        url += "?" + urlencode(params)
    request = Request(url, headers={
        "Accept": "application/json",
        "User-Agent": "MVA-WGS-reference-context/1.0",
    })
    last_error = None
    for attempt in range(3):
        try:
            with urlopen(request, timeout=90) as response:
                raw = response.read()
            return json.loads(raw), raw, url
        except Exception as error:
            last_error = error
            if attempt < 2:
                time.sleep(2 ** (attempt + 1))
    raise RuntimeError(f"Reference request failed: {url}") from last_error

def load_reference(server, transcript, version, gene, cache_root, allow_live):
    host_key = sha(server.encode())[:12]
    directory = (
        Path(cache_root) / host_key / f"{transcript}.{version}"
    )
    manifest_path = directory / "manifest.json"

    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        for filename, expected in manifest["sha256"].items():
            if sha((directory / filename).read_bytes()) != expected:
                raise ValueError("Reference-cache checksum mismatch")
        meta = json.loads((directory / "transcript.json").read_bytes())
        sequence = json.loads((directory / "cds.json").read_bytes())
        gene_meta = json.loads((directory / "gene.json").read_bytes())
    else:
        if directory.exists():
            raise ValueError("Incomplete reference cache; inspect before reuse")
        if not allow_live:
            raise ValueError("Reference not cached and live retrieval disabled")

        meta, meta_raw, meta_url = fetch_json(
            server, "/lookup/id/" + quote(transcript, safe="")
        )
        validate_lookup(meta, transcript, version)

        parent = meta.get("Parent")
        if not parent or not re.fullmatch(r"ENSG[0-9]+", parent):
            raise ValueError("Transcript lacks a valid parent gene ID")

        gene_meta, gene_raw, gene_url = fetch_json(
            server, "/lookup/id/" + quote(parent, safe="")
        )
        if gene_meta.get("display_name") != gene:
            raise ValueError("Reference transcript belongs to another gene")

        sequence, sequence_raw, sequence_url = fetch_json(
            server,
            "/sequence/id/" + quote(transcript, safe=""),
            {"type": "cds"},
        )
        seq_id = str(sequence.get("id", ""))
        if seq_id not in {transcript, f"{transcript}.{version}"}:
            raise ValueError("Unexpected CDS identifier")

        # Check the live transcript version again after sequence acquisition.
        after, after_raw, after_url = fetch_json(
            server, "/lookup/id/" + quote(transcript, safe="")
        )
        validate_lookup(after, transcript, version)
        if after.get("Parent") != parent:
            raise ValueError("Transcript parent changed during retrieval")

        payloads = {
            "transcript.json": meta_raw,
            "gene.json": gene_raw,
            "cds.json": sequence_raw,
            "transcript_after.json": after_raw,
        }
        manifest = {
            "reference_id": f"{transcript}.{version}",
            "gene_symbol": gene,
            "server": server,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "assembly_name": meta.get("assembly_name"),
            "numeric_ensembl_release": None,
            "urls": [meta_url, gene_url, sequence_url, after_url],
            "sha256": {
                filename: sha(raw) for filename, raw in payloads.items()
            },
        }
        directory.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=".reference_stage_", dir=directory.parent
        ) as temporary:
            stage = Path(temporary)
            for filename, raw in payloads.items():
                (stage / filename).write_bytes(raw)
            (stage / "manifest.json").write_text(
                json.dumps(manifest, indent=2) + "\n"
            )
            os.rename(stage, directory)

    validate_lookup(meta, transcript, version)
    if gene_meta.get("display_name") != gene:
        raise ValueError("Cached transcript/gene mismatch")
    seq_id = str(sequence.get("id", ""))
    if seq_id not in {transcript, f"{transcript}.{version}"}:
        raise ValueError("Cached CDS identifier mismatch")
    cds = str(sequence.get("seq", "")).upper()
    if not cds or set(cds) - set("ACGTN"):
        raise ValueError("Invalid reference CDS")

    return cds, meta, directory

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--a10", required=True)
    parser.add_argument("--gene", action="append", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--server", default="https://rest.ensembl.org")
    parser.add_argument("--flank", type=int, default=15)
    parser.add_argument("--allow-live", action="store_true")
    args = parser.parse_args()

    if args.flank < 1:
        raise ValueError("Flank length must be positive")
    if not args.server.startswith("https://"):
        raise ValueError("HTTPS reference server required")

    evidence = read_tsv(args.evidence)
    a10 = read_tsv(args.a10)
    selected_genes = set(args.gene)

    representatives = {}
    for row in a10:
        gene = row.get("gene_symbol", "")
        if gene in selected_genes:
            if gene in representatives:
                raise ValueError("Duplicate A10 representative gene")
            representatives[gene] = row

    if set(representatives) != selected_genes:
        raise ValueError("Configured gene missing from A10")

    selected = [
        row for row in evidence
        if row.get("gene_symbol", row.get("gene", "")) in selected_genes
    ]
    if not selected:
        raise ValueError("No retained variant evidence for configured genes")

    references = {}
    output = []
    seen = set()

    for row in selected:
        gene = row.get("gene_symbol", row.get("gene", ""))
        hgvsc, hgvsp = row.get("HGVSc", ""), row.get("HGVSp", "")
        v = parse_hgvs(hgvsc, hgvsp)

        feature = row.get("Feature", "").split(".")[0]
        if feature and feature != v["transcript"]:
            raise ValueError("Feature and versioned HGVSc disagree")

        candidate = row["candidate_id"]
        key = (gene, candidate, hgvsc, hgvsp)
        if key in seen:
            continue
        seen.add(key)

        representative = (
            candidate == representatives[gene]["best_candidate_id"]
        )
        if representative:
            a10_tx = representatives[gene].get("transcript_id", "").split(".")[0]
            if a10_tx and a10_tx != v["transcript"]:
                raise ValueError("A10 and retained representative transcripts differ")

        ref_key = (v["transcript"], v["version"], gene)
        if ref_key not in references:
            references[ref_key] = load_reference(
                args.server, v["transcript"], v["version"], gene,
                args.cache, args.allow_live,
            )
        cds, meta, reference_dir = references[ref_key]
        result = annotate(cds, hgvsc, hgvsp, args.flank)

        gt = row.get("GT", "")
        result = {
            "gene": gene,
            "candidate_id": candidate,
            "GT": gt,
            "phase_status": (
                "unphased" if "/" in gt else
                "phased_GT_present_pair_phase_not_validated" if "|" in gt else
                "not_established"
            ),
            "consequence": row.get("consequence", ""),
            "representative_in_A10": representative,
            **result,
            "reference_assembly": meta.get("assembly_name", ""),
            "reference_manifest": str(reference_dir / "manifest.json"),
            "reference_CDS_sha256": sha(cds.encode()),
        }
        output.append(result)

    for gene in selected_genes:
        if not any(
            row["gene"] == gene and row["representative_in_A10"]
            for row in output
        ):
            raise ValueError("A10 representative absent from context output")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    table = out / "codon_context.tsv"
    fasta = out / "flanking_context.fasta"

    with table.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(output[0]),
            delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(output)

    with fasta.open("w", encoding="utf-8") as stream:
        for i, row in enumerate(output, 1):
            for label, field in (
                ("reference_DNA", "reference_context_DNA"),
                ("focal_variant_DNA", "alternate_context_DNA"),
                ("reference_RNA", "reference_context_RNA"),
                ("focal_variant_RNA", "alternate_context_RNA"),
            ):
                stream.write(
                    f">context_{i}|{row['gene']}|"
                    f"{row['transcript_version']}|{label}\n"
                    f"{row[field]}\n"
                )

    manifest = {
        "node": "d3_codon_context",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": sha(Path(__file__).read_bytes()),
        "input_sha256": {
            args.evidence: sha(Path(args.evidence).read_bytes()),
            args.a10: sha(Path(args.a10).read_bytes()),
        },
        "selected_genes": sorted(selected_genes),
        "requested_flank_bases_each_side": args.flank,
        "reference_manifests": {
            str(directory / "manifest.json"):
                sha((directory / "manifest.json").read_bytes())
            for _, _, directory in references.values()
        },
        "output_sha256": {
            table.name: sha(table.read_bytes()),
            fasta.name: sha(fasta.read_bytes()),
        },
        "sequence_verified": True,
        "genomic_coordinate_mapping_verified": False,
        "phase_established": False,
        "NMD_assessed": False,
        "compound_efficacy_assessed": False,
        "remote_requests_contain": "public_reference_identifiers_only",
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    print("CODON CONTEXT COMPLETE")
    for row in output:
        print(
            row["gene"], row["HGVSc"], row["HGVSp"],
            row["reference_codon_DNA"] + ">" + row["alternate_codon_DNA"],
            row["alternate_context_bracketed_RNA"],
            "stop_plus4=" + row["stop_plus4_motif_RNA"],
            sep="\t",
        )

if __name__ == "__main__":
    main()
