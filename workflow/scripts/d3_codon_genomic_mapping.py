#!/usr/bin/env python3
"""Fail-closed, positive-strand GRCh38 CDS-context mapping producer."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import tempfile
import time
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

IDENTITY = ["gene", "candidate_id", "transcript_version", "HGVSc", "HGVSp",
            "alternate_codon_RNA", "alternate_context_bracketed_RNA"]
FIELDS = IDENTITY + ["assembly", "chrom", "strand", "bed_chromStart", "bed_chromEnd",
    "genomic_coordinate_mapping_verified", "genomic_context_matches_frozen_CDS",
    "codon_contiguous", "track_lookup_status", "cds_context_start", "cds_context_end",
    "genomic_reference_context_DNA", "validation_scope"]

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def require(condition, message):
    if not condition:
        raise ValueError(message)

def verify_outputs(directory, manifest):
    hashes = manifest.get("output_sha256")
    require(isinstance(hashes, dict) and bool(hashes), "Missing output hashes")
    tracked = {}
    for name, expected in hashes.items():
        relative = Path(name)
        require(not relative.is_absolute() and ".." not in relative.parts,
                "Unsafe manifest path")
        path = directory / relative
        require(sha(path.read_bytes()) == expected, f"Checksum mismatch: {path}")
        tracked[str(path)] = expected
    return tracked

def verify_transcript(meta, row):
    match = re.fullmatch(r"(ENST\d+)\.(\d+)", row["transcript_version"])
    require(match is not None, "Expected versioned ENST identifier")
    tx, version = match.groups()
    require(meta.get("id") == tx and int(meta.get("version", -1)) == int(version),
            "Reference transcript version mismatch; use a matching archive")
    require(meta.get("species") == "homo_sapiens" and
            meta.get("object_type") == "Transcript" and
            meta.get("assembly_name") == "GRCh38", "Transcript reference scope mismatch")
    return tx

def map_record(row, response, sequence):
    require(str(row.get("sequence_verified")).lower() == "true", "Unverified CDS context")
    require(str(row.get("alternate_is_stop")).lower() == "true", "Not a stop record")
    require(str(row.get("context_boundary_truncated")).lower() == "false", "Truncated context")
    require(row.get("reference_assembly") == "GRCh38", "Context assembly mismatch")
    h = re.fullmatch(r"(ENST\d+\.\d+):c\.(\d+)([ACGT])>([ACGT])", row["HGVSc"])
    require(h is not None, "Only versioned coding SNVs supported")
    tx, position, ref, alt = h.groups(); position = int(position)
    require(tx == row["transcript_version"], "HGVSc transcript mismatch")
    codon_start = ((position - 1) // 3) * 3 + 1
    left, right = int(row["upstream_bases"]), int(row["downstream_bases"])
    require(left > 0 and right > 0, "Missing flanks")
    start, end = codon_start - left, codon_start + 2 + right
    require(start >= 1, "Invalid CDS interval")
    reference = row["reference_context_DNA"].upper()
    require(len(reference) == end-start+1 and not set(reference)-set("ACGT"),
            "Malformed reference window")
    require(reference[position-start] == ref, "HGVSc REF differs from CDS")
    mutant = reference[:position-start] + alt + reference[position-start+1:]
    require(mutant == row["alternate_context_DNA"], "Mutant window disagreement")
    stop = mutant[left:left+3].replace("T", "U")
    require(stop in {"UGA", "UAG", "UAA"} and stop == row["alternate_codon_RNA"],
            "Altered codon is not the recorded stop")
    bracket = mutant[:left].replace("T", "U") + "[" + stop + "]" + mutant[left+3:].replace("T", "U")
    require(bracket == row["alternate_context_bracketed_RNA"], "RNA context disagreement")
    blocks = response.get("mappings")
    require(isinstance(blocks, list) and bool(blocks), "No mapping blocks")
    coordinates, bases, chromosomes = {}, {}, set()
    for block in blocks:
        original, mapped = block["original"], block["mapped"]
        a, b = int(original["start"]), int(original["end"])
        gs, ge = int(mapped["start"]), int(mapped["end"])
        require(start <= a <= b <= end, "Mapping block outside requested CDS window")
        require(mapped.get("assembly_name", mapped.get("assembly")) == "GRCh38",
                "Mapping assembly mismatch")
        require(int(mapped["strand"]) == 1, "Negative strand unsupported by current downstream adapter")
        require(ge-gs == b-a, "Gapped or inconsistent mapping block")
        chrom = str(mapped["seq_region_name"])
        require(re.fullmatch(r"(?:[1-9]|1\d|2[0-2]|X|Y|MT)", chrom) is not None,
                "Unsupported chromosome/contig")
        chrom = "chrM" if chrom == "MT" else "chr" + chrom
        chromosomes.add(chrom)
        dna = sequence(chrom, gs-1, ge).upper()
        require(len(dna) == b-a+1 and not set(dna)-set("ACGT"), "Invalid genomic sequence block")
        for offset, cds_position in enumerate(range(a,b+1)):
            require(cds_position not in coordinates, "Overlapping/ambiguous CDS mappings")
            coordinates[cds_position] = gs + offset
            bases[cds_position] = dna[offset]
    require(len(chromosomes) == 1, "Multiple genomic chromosomes")
    require(set(coordinates) == set(range(start,end+1)), "Incomplete CDS-window mapping")
    genomic = "".join(bases[p] for p in range(start,end+1))
    require(genomic == reference, "Full genomic reference context differs from frozen CDS")
    variant = row["candidate_id"].split(":")
    require(len(variant) == 4, "Expected CHROM:POS:REF:ALT candidate key")
    chrom, genomic_position, genomic_ref, genomic_alt = variant
    expected_chrom = next(iter(chromosomes))
    require(chrom.removeprefix("chr") == expected_chrom.removeprefix("chr"), "Candidate chromosome mismatch")
    require(int(genomic_position) == coordinates[position], "Candidate genomic position mismatch")
    require(genomic_ref == ref and genomic_alt == alt, "Candidate alleles disagree with coding SNV")
    codon = [coordinates[p] for p in range(codon_start,codon_start+3)]
    require(codon == list(range(codon[0],codon[0]+3)), "Stop codon spans a genomic discontinuity")
    result = {key: row[key] for key in IDENTITY}
    result.update(assembly="GRCh38", chrom=expected_chrom, strand="1",
        bed_chromStart=str(codon[0]-1), bed_chromEnd=str(codon[-1]),
        genomic_coordinate_mapping_verified="True", genomic_context_matches_frozen_CDS="True",
        codon_contiguous="True", track_lookup_status="contiguous_stop_codon_lookup_key",
        cds_context_start=str(start), cds_context_end=str(end),
        genomic_reference_context_DNA=genomic,
        validation_scope="full_reference_CDS_window_and_exact_candidate_SNV_positive_strand")
    return result

class Requests:
    def __init__(self, directory, live):
        self.directory, self.live, self.used = directory, live, {}
        self.network_requests = 0
    def get(self, url):
        path = self.directory / (sha(url.encode()) + ".json")
        if path.exists():
            envelope = json.loads(path.read_text())
            require(envelope["url"] == url, "Cached URL mismatch")
            raw = envelope["raw_response"].encode()
            require(sha(raw) == envelope["response_sha256"], "Request-cache checksum mismatch")
        else:
            require(self.live, f"Request not cached; live retrieval disabled: {url}")
            raw = None
            for attempt in range(3):
                try:
                    with urlopen(Request(url, headers={"Accept":"application/json",
                        "User-Agent":"MVA-CDS-context-mapping/1"}), timeout=90) as response:
                        raw = response.read()
                    self.network_requests += 1
                    break
                except Exception:
                    if attempt == 2:
                        raise
                    time.sleep(2**attempt)
            json.loads(raw)
            envelope = {"url":url,"retrieved_utc":datetime.now(timezone.utc).isoformat(),
                        "response_sha256":sha(raw),"raw_response":raw.decode("utf-8")}
            self.directory.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(mode="w", dir=self.directory, delete=False) as stream:
                json.dump(envelope, stream); temporary = Path(stream.name)
            temporary.replace(path)
        self.used[str(path)] = sha(path.read_bytes())
        return json.loads(raw)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--request-cache", type=Path, required=True)
    parser.add_argument("--server", default="https://rest.ensembl.org")
    parser.add_argument("--allow-live", action="store_true")
    args = parser.parse_args()
    require(args.server.startswith("https://"), "HTTPS server required")
    require(args.outdir.resolve() != args.context_dir.resolve(), "Separate output directory required")
    context = args.context_dir.resolve()
    manifest_path = context / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    tracked = {str(manifest_path):sha(manifest_path.read_bytes())}
    tracked.update(verify_outputs(context, manifest))
    require("codon_context.tsv" in manifest["output_sha256"], "Context table not pinned")
    for key in ("input_sha256", "reference_manifests"):
        entries = manifest.get(key)
        require(isinstance(entries, dict) and bool(entries), f"Missing provenance: {key}")
        for filename, expected in entries.items():
            require(sha(Path(filename).read_bytes()) == expected, f"Upstream fingerprint mismatch: {filename}")
            tracked[filename] = expected
    with (context / "codon_context.tsv").open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames and len(reader.fieldnames)==len(set(reader.fieldnames)), "Invalid table headers")
        rows = list(reader)
    require(all(None not in row and None not in row.values() for row in rows), "Malformed context rows")
    selected = [row for row in rows if str(row.get("alternate_is_stop")).lower()=="true"]
    require(len({(r["gene"],r["candidate_id"]) for r in selected})==len(selected), "Ambiguous stop contexts")
    requests = Requests(args.request_cache.resolve(), args.allow_live)
    results = []
    for row in selected:
        reference_path = Path(row["reference_manifest"])
        require(str(reference_path) in manifest["reference_manifests"], "Unpinned per-record reference")
        reference_manifest = json.loads(reference_path.read_text())
        for name, expected in reference_manifest["sha256"].items():
            rel = Path(name)
            require(not rel.is_absolute() and ".." not in rel.parts, "Unsafe reference path")
            path = reference_path.parent / rel
            require(sha(path.read_bytes()) == expected, "Frozen reference checksum mismatch")
            tracked[str(path)] = expected
        frozen = json.loads((reference_path.parent / "transcript.json").read_text())
        tx = verify_transcript(frozen, row)
        cds = json.loads((reference_path.parent / "cds.json").read_text())["seq"].upper()
        require(sha(cds.encode()) == row["reference_CDS_sha256"], "Frozen CDS fingerprint mismatch")
        coding_position = int(row["coding_position"])
        parsed_position = re.fullmatch(r"ENST\d+\.\d+:c\.(\d+)[ACGT]>[ACGT]", row["HGVSc"])
        require(parsed_position is not None and int(parsed_position.group(1)) == coding_position,
                "Coding position disagrees with HGVSc")
        cstart = ((coding_position-1)//3)*3+1
        begin = cstart-int(row["upstream_bases"])
        finish = cstart+2+int(row["downstream_bases"])
        require(cds[begin-1:finish] == row["reference_context_DNA"], "Context is not the pinned CDS slice")
        base = args.server.rstrip("/")
        lookup_url = base + "/lookup/id/" + quote(tx, safe="")
        before = requests.get(lookup_url)
        verify_transcript(before, row)
        require(before.get("Parent") == frozen.get("Parent"), "Transcript gene changed")
        mapping = requests.get(base + "/map/cds/" + tx + f"/{begin}..{finish}?include_original_region=1")
        def sequence(chrom, start, end):
            url = "https://api.genome.ucsc.edu/getData/sequence?" + urlencode(
                {"genome":"hg38","chrom":chrom,"start":start,"end":end})
            payload = requests.get(url)
            require(payload.get("genome")=="hg38" and payload.get("chrom")==chrom,
                    "Genomic sequence response identity mismatch")
            require(int(payload["start"])==start and int(payload["end"])==end,
                    "Genomic sequence interval mismatch")
            return payload["dna"]
        results.append(map_record(row,mapping,sequence))
        after = requests.get(lookup_url + "?expand=0")
        verify_transcript(after,row)
        require(after.get("Parent")==before.get("Parent"), "Transcript identity changed during retrieval")
    tracked.update(requests.used)
    tracked[str(Path(__file__).resolve())] = sha(Path(__file__).read_bytes())
    for filename, expected in tracked.items():
        require(sha(Path(filename).read_bytes())==expected, f"Input changed: {filename}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".mapping_", dir=args.outdir) as temporary:
        stage = Path(temporary)
        with (stage / "codon_genomic_mapping.tsv").open("w",newline="") as stream:
            writer = csv.DictWriter(stream,fieldnames=FIELDS,delimiter="\t",lineterminator="\n")
            writer.writeheader(); writer.writerows(results)
        with (stage / "stop_codon_lookup.bed").open("w") as stream:
            for r in results:
                stream.write("\t".join([r["chrom"],r["bed_chromStart"],r["bed_chromEnd"],r["candidate_id"],"0","+"])+"\n")
        audit = {"node":"stop_codon_genomic_mapping_audit_v2",
            "created_utc":datetime.now(timezone.utc).isoformat(),
            "source_context_records":len(rows),"selected_stop_records":len(selected),
            "mapping_scope":"full_CDS_context_validation_then_three_base_stop_lookup",
            "supported_scope":"GRCh38_positive_strand_coding_SNV_contiguous_stop_codon",
            "input_and_implementation_sha256":tracked,
            "network_allowed":args.allow_live,"network_requests_made":requests.network_requests,
            "request_cache_files":requests.used,
            "transcript_version_checked_in_mapping_server_responses":True,
            "mapping_endpoint_version_binding_independently_guaranteed":False,
            "clinical_eligibility_assessed":False,
            "output_sha256":{name:sha((stage/name).read_bytes()) for name in
                ["codon_genomic_mapping.tsv","stop_codon_lookup.bed"]}}
        (stage/"manifest.json").write_text(json.dumps(audit,indent=2)+"\n")
        for name in ("codon_genomic_mapping.tsv","stop_codon_lookup.bed","manifest.json"):
            (stage/name).replace(args.outdir/name)
    print(json.dumps({"mapped_stop_records":len(results),"outdir":str(args.outdir)}))

if __name__ == "__main__":
    main()
