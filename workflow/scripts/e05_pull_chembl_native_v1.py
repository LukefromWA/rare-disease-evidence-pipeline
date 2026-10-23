#!/usr/bin/env python3
"""Acquire ChEMBL mechanisms and referenced targets/molecules, without biological normalization."""
import argparse
import collections
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import urllib.parse as up

BASE = "https://www.ebi.ac.uk/chembl/api/data/"
LICENSE = "CC BY-SA 3.0 Unported"

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def write_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)

def api_url(endpoint, **params):
    return BASE + endpoint + ".json?" + up.urlencode(params)

def checked_url(url):
    url = up.urljoin(BASE, url)
    p = up.urlsplit(url)
    if p.scheme != "https" or p.netloc != "www.ebi.ac.uk" or not p.path.startswith("/chembl/api/data/"):
        raise RuntimeError("Unexpected API pagination URL: " + url)
    return url

class Client:
    def __init__(self, out, curl, delay):
        self.out, self.curl, self.delay = Path(out), curl, delay

    def get(self, url, label):
        url = checked_url(url)
        path = self.out / "pages" / (label + ".json")
        meta = path.with_suffix(".request.json")
        if path.exists() or meta.exists():
            if not (path.is_file() and meta.is_file()):
                raise RuntimeError("Incomplete cached request; inspect before removing: " + str(path))
            saved = json.loads(meta.read_text())
            if saved["url"] != url or saved["sha256"] != sha(path):
                raise RuntimeError("Cached URL/hash mismatch: " + str(path))
            print("CACHE", label, flush=True)
            return json.loads(path.read_text())
        part = path.with_suffix(".json.part")
        common = [self.curl, "--proto", "=https", "--proto-redir", "=https", "--http1.1",
                  "--fail", "--silent", "--show-error", "--location", "--connect-timeout", "20",
                  "--max-time", "120", "--retry", "2", "--retry-delay", "2", "--retry-all-errors"]
        for mode, extra in [("default", []), ("tls12", ["--tlsv1.2", "--tls-max", "1.2"])]:
            headers = path.with_name(label + "." + mode + ".headers.txt")
            err = path.with_name(label + "." + mode + ".stderr.txt")
            command = common + extra + ["--dump-header", str(headers), "--output", str(part), url]
            print("GET", label, mode, flush=True)
            with err.open("wb") as stream:
                result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=stream)
            if result.returncode == 0:
                obj = json.loads(part.read_text(encoding="utf-8"))
                if not isinstance(obj, dict):
                    raise RuntimeError("Expected JSON object: " + label)
                os.replace(part, path)
                write_json(meta, {"url": url, "retrieved_utc": utc(), "transport": mode,
                                 "sha256": sha(path), "headers_sha256": sha(headers)})
                time.sleep(self.delay)
                return obj
            print(err.read_text(errors="replace").strip(), file=sys.stderr, flush=True)
        raise RuntimeError("Both verified-HTTPS transport attempts failed: " + label)

def crawl(client, url, key, id_key, prefix):
    records, ids, visited = [], set(), set()
    total = None
    page = 0
    while url:
        url = checked_url(url)
        if url in visited:
            raise RuntimeError("Pagination cycle: " + url)
        visited.add(url)
        obj = client.get(url, f"{prefix}_{page:04d}")
        batch, meta = obj.get(key), obj.get("page_meta")
        if not isinstance(batch, list) or not isinstance(meta, dict):
            raise RuntimeError("Unexpected collection schema: " + prefix)
        count = meta.get("total_count")
        if not isinstance(count, int) or isinstance(count, bool) or count < 0:
            raise RuntimeError("Missing/invalid collection count")
        if total is None:
            total = count
        if count != total or meta.get("offset") != len(records):
            raise RuntimeError("Collection count/offset changed: " + prefix)
        for rec in batch:
            identity = rec.get(id_key)
            if identity is None or identity in ids:
                raise RuntimeError("Missing/duplicate collection identity: " + prefix)
            ids.add(identity)
            records.append(rec)
        print(prefix, f"{len(records)}/{total}", flush=True)
        nxt = meta.get("next")
        if nxt and not batch:
            raise RuntimeError("Empty nonfinal page")
        url = checked_url(nxt) if nxt else None
        page += 1
    if len(records) != total:
        raise RuntimeError("Incomplete collection: " + prefix)
    return records

def referenced(client, endpoint, key, id_key, ids, batch_size, fields=None):
    found = {}
    ids = sorted(ids)
    for index in range(0, len(ids), batch_size):
        wanted = set(ids[index:index+batch_size])
        query = {id_key + "__in": ",".join(sorted(wanted)), "limit": 1000, "offset": 0, "order_by": id_key}
        if fields:
            query["only"] = ",".join(fields)
        rows = crawl(client, api_url(endpoint, **query), key, id_key, f"{endpoint}_batch{index//batch_size:04d}")
        got = {r[id_key] for r in rows}
        if got != wanted:
            raise RuntimeError(f"{endpoint} lookup mismatch; missing={sorted(wanted-got)} unexpected={sorted(got-wanted)}")
        if set(found).intersection(got):
            raise RuntimeError("Duplicate entity across batches")
        found.update({r[id_key]: r for r in rows})
    return list(found.values())

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True)
    p.add_argument("--expected-version", default="ChEMBL_37")
    p.add_argument("--expected-mechanisms", type=int, default=7561)
    p.add_argument("--batch-size", type=int, default=50)
    p.add_argument("--delay", type=float, default=0.25)
    p.add_argument("--curl", default="/usr/bin/curl")
    p.add_argument("--resume", action="store_true")
    a = p.parse_args()
    if not 1 <= a.batch_size <= 100 or a.delay < 0:
        p.error("batch size must be 1..100 and delay must be nonnegative")
    curl = shutil.which(a.curl)
    if not curl:
        p.error("curl executable not found")
    out = Path(a.out).resolve()
    if (out / "download_manifest.json").exists():
        p.error("Acquisition already complete; verify native.sha256 rather than overwrite")
    if a.resume:
        if not out.is_dir():
            p.error("--resume requires the exact existing incomplete output directory")
    else:
        out.mkdir(parents=True, exist_ok=False)
    (out / "pages").mkdir(exist_ok=True)
    client = Client(out, curl, a.delay)
    token = str(time.time_ns())
    # Always re-query status on resume; cached resource pages require exact URL/hash matches.
    before = client.get(BASE + "status.json", "status_start_" + token)
    if before.get("chembl_db_version") != a.expected_version or before.get("status") != "UP":
        raise RuntimeError("Unexpected ChEMBL version or API status: " + json.dumps(before))
    mechanism = crawl(client, api_url("mechanism", limit=1000, offset=0, order_by="mec_id"),
                      "mechanisms", "mec_id", "mechanism")
    if len(mechanism) != a.expected_mechanisms:
        raise RuntimeError(f"Expected {a.expected_mechanisms} mechanisms, received {len(mechanism)}; review before accepting")
    target_ids = {r["target_chembl_id"] for r in mechanism if r.get("target_chembl_id")}
    molecule_ids = {r[k] for r in mechanism for k in ("molecule_chembl_id", "parent_molecule_chembl_id") if r.get(k)}
    if any(not re.fullmatch(r"CHEMBL[0-9]+", x) for x in target_ids | molecule_ids):
        raise RuntimeError("Unexpected ChEMBL identifier format")
    targets = referenced(client, "target", "targets", "target_chembl_id", target_ids, a.batch_size,
                         ["target_chembl_id", "pref_name", "target_type", "organism", "tax_id", "species_group_flag", "target_components"])
    molecules = referenced(client, "molecule", "molecules", "molecule_chembl_id", molecule_ids, a.batch_size,
                           ["molecule_chembl_id", "pref_name", "molecule_type", "max_phase"])
    after = client.get(BASE + "status.json", "status_end_" + token)
    if after.get("chembl_db_version") != a.expected_version or after.get("chembl_release_date") != before.get("chembl_release_date") or after.get("status") != "UP":
        raise RuntimeError("API version/release/status changed during acquisition")
    action_counts = collections.Counter(str(r.get("action_type")) for r in mechanism)
    type_counts = collections.Counter(str(r.get("target_type")) for r in targets)
    manifest = {"format_version": 1, "finished_utc": utc(), "source": "ChEMBL", "api_base": BASE,
        "chembl_db_version": a.expected_version, "chembl_release_date": before.get("chembl_release_date"),
        "license": LICENSE, "license_documentation": "https://chembl.gitbook.io/chembl-interface-documentation/about",
        "counts": {"mechanisms": len(mechanism), "referenced_targets": len(targets), "referenced_molecules_and_parents": len(molecules),
                   "mechanisms_without_target": sum(not r.get("target_chembl_id") for r in mechanism)},
        "action_types": dict(sorted(action_counts.items())), "target_types": dict(sorted(type_counts.items())),
        "molecule_fields_requested": ["molecule_chembl_id", "pref_name", "molecule_type", "max_phase"],
        "scope": "All mechanism records, referenced targets with nested components, referenced molecules and parents. No approval, human, clinical-efficacy, direct-interaction or target-type filtering applied during acquisition.",
        "limitations": "Live API version/count checks are not a transactionally frozen database snapshot. No signed gene edges, gene mapping, model fitting, clinical ranking or readiness certification performed.",
        "downloader_sha256": sha(__file__), "downloader_filename": Path(__file__).name,
        "curl_version": subprocess.check_output([curl, "--version"], text=True).splitlines()[0]}
    code_dir = out / "acquisition_code"
    code_dir.mkdir(exist_ok=True)
    shutil.copyfile(__file__, code_dir / Path(__file__).name)
    # Write checksums before the completion manifest; existence of the manifest means these checks finished.
    contents = [x for x in out.rglob("*") if x.is_file() and x.name not in ("native.sha256", "download_manifest.json") and not x.name.endswith((".part", ".tmp"))]
    checksum_text = "".join(sha(x) + "  " + x.relative_to(out).as_posix() + "\n" for x in sorted(contents))
    (out / "native.sha256").write_text(checksum_text)
    manifest["native_sha256_file_sha256"] = sha(out / "native.sha256")
    write_json(out / "download_manifest.json", manifest)
    print("\nCHEMBL NATIVE ACQUISITION COMPLETE")
    print(json.dumps(manifest["counts"], indent=2))
    print("\nACTION TYPES:", json.dumps(manifest["action_types"], indent=2))
    print("\nTARGET TYPES:", json.dumps(manifest["target_types"], indent=2))
    print("\nOutput:", out)
    print("No normalized drug-target table, source/config edits, fitting or readiness flags.")

if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, OSError, subprocess.SubprocessError) as exc:
        print("STOP:", str(exc), file=sys.stderr)
        print("Retain the incomplete directory. Resume with the same arguments plus --resume; do not change batching parameters.", file=sys.stderr)
        sys.exit(1)
