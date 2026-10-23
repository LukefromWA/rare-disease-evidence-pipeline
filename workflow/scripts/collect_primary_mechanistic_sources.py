#!/usr/bin/env python3
import argparse
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import re
import tempfile
import time
import unicodedata
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

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
PASSAGE_FIELDS = [
    "source_id", "pmcid", "section_locator",
    "paragraph_number", "matched_topics", "source_text", "review_status",
]

TOPICS = {
    "termination_mechanism": [
        r"release.factor", r"near.cognate", r"termination",
        r"purified", r"kinetic",
    ],
    "model_context": [
        r"upstream", r"downstream", r"interaction term",
        r"logistic", r"sequence context", r"prediction model",
    ],
    "endpoint_and_normalization": [
        r"normaliz", r"percentage", r"percent\b", r"no.nonsense",
        r"reporter", r"fluorescen",
    ],
    "functional_limitations": [
        r"not necessarily", r"functionality", r"functional protein",
        r"amino.acid", r"substitut", r"rescue",
    ],
}

def sha(data):
    return hashlib.sha256(data).hexdigest()

def norm(value):
    value = unicodedata.normalize("NFKC", str(value))
    for character in ["−", "–", "—"]:
        value = value.replace(character, "-")
    return " ".join(value.split())

def text(element):
    if element is None:
        return ""
    return " ".join("".join(element.itertext()).split())

def write_json(path, value):
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fields,
            delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

def parse_article(raw, expected):
    root = ET.fromstring(raw)
    candidates = (
        [root] if root.tag == "article"
        else root.findall(".//article")
    )
    matches = []

    for article in candidates:
        meta = article.find("./front/article-meta")
        if meta is None:
            continue
        identifiers = {}
        for element in meta.findall("article-id"):
            identifiers.setdefault(
                element.get("pub-id-type", ""), []
            ).append(text(element))

        pmc_ids = {
            value.upper().removeprefix("PMC")
            for key in ["pmc", "pmcid"]
            for value in identifiers.get(key, [])
        }
        dois = {
            value.lower() for value in identifiers.get("doi", [])
        }
        pmids = set(identifiers.get("pmid", []))

        if (
            expected["pmcid"].upper().removeprefix("PMC") not in pmc_ids
            or expected["doi"].lower() not in dois
        ):
            continue
        if expected.get("pmid") and expected["pmid"] not in pmids:
            continue
        matches.append((article, meta, identifiers))

    if len(matches) != 1:
        raise ValueError("PMC/DOI/PMID identity was not uniquely verified")

    article, meta, identifiers = matches[0]
    body = article.find("body")
    if body is None or not text(body):
        raise ValueError("Response has no usable full-text body")

    pmids = identifiers.get("pmid", [])
    if len(set(pmids)) != 1:
        raise ValueError("Article lacks a unique PMID")
    pmid = pmids[0]

    authors = []
    bib_authors = []
    for contributor in meta.findall(".//contrib"):
        if contributor.get("contrib-type") != "author":
            continue
        name = contributor.find("name")
        if name is not None:
            surname = text(name.find("surname"))
            given = text(name.find("given-names"))
            if surname:
                authors.append(f"{surname}, {given}".strip(", "))
                bib_authors.append(f"{surname}, {given}".strip(", "))
        else:
            collective = text(contributor.find("collab"))
            if collective:
                authors.append(collective)
                bib_authors.append("{" + collective + "}")

    year = ""
    for date in meta.findall("pub-date"):
        value = text(date.find("year"))
        if value:
            year = value
            break

    abstract_parts = [
        text(element) for element in meta.findall("abstract")
        if text(element)
    ]

    info = {
        "source_id": "PMID" + pmid,
        "pmid": pmid,
        "pmcid": expected["pmcid"],
        "doi": expected["doi"],
        "role": expected["role"],
        "title": text(meta.find("title-group/article-title")),
        "authors": "; ".join(authors),
        "bib_authors": bib_authors,
        "year": year,
        "publication": text(article.find(
            "front/journal-meta/journal-title-group/journal-title"
        )),
        "volume": text(meta.find("volume")),
        "issue": text(meta.find("issue")),
        "pages_or_article_id": (
            text(meta.find("elocation-id"))
            or text(meta.find("fpage"))
        ),
        "abstract": "\n\n".join(abstract_parts),
        "license_text": text(meta.find("permissions/license")),
    }
    return article, info

def collect_passages(article, info):
    parents = {
        child: parent
        for parent in article.iter()
        for child in parent
    }

    def location(element):
        parts = []
        current = element
        while current is not None:
            if current.tag == "ref-list":
                return None
            if current.tag == "sec":
                parts.append(
                    (
                        current.get("id", "") + ": "
                        + text(current.find("title"))
                    ).strip(": ")
                )
            elif current.tag in {"abstract", "body", "back"}:
                parts.append(current.tag)
            current = parents.get(current)
        return " > ".join(reversed(parts))

    rows = []
    paragraph_number = 0
    for paragraph in article.iter("p"):
        locator = location(paragraph)
        if locator is None:
            continue
        value = text(paragraph)
        if not value:
            continue
        paragraph_number += 1
        topics = [
            topic for topic, patterns in TOPICS.items()
            if any(re.search(pattern, value, re.I) for pattern in patterns)
        ]
        if topics:
            rows.append({
                "source_id": info["source_id"],
                "pmcid": info["pmcid"],
                "section_locator": locator,
                "paragraph_number": paragraph_number,
                "matched_topics": ";".join(topics),
                "source_text": value,
                "review_status": "passage_selected_not_adjudicated",
            })
    return rows

def cache_article(expected, cache_root):
    directory = cache_root / expected["pmcid"]

    if directory.exists():
        manifest_path = directory / "manifest.json"
        xml_path = directory / "article.xml"
        if not manifest_path.is_file() or not xml_path.is_file():
            raise ValueError(
                f"Incomplete cache directory; inspect manually: {directory}"
            )
        manifest = json.loads(manifest_path.read_text())
        raw = xml_path.read_bytes()
        if (
            manifest.get("expected_identity") != expected
            or manifest.get("article_sha256") != sha(raw)
        ):
            raise ValueError("Existing article cache provenance mismatch")
        article, info = parse_article(raw, expected)
        print("Reusing verified:", expected["pmcid"], flush=True)
        return raw, article, info, manifest

    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        + urlencode({
            "db": "pmc",
            "id": expected["pmcid"].removeprefix("PMC"),
            "retmode": "xml",
            "tool": "primary_mechanistic_source_collection",
        })
    )
    last_error = None
    for attempt in range(3):
        time.sleep(0.5)
        try:
            print("Fetching:", expected["pmcid"], flush=True)
            request = Request(
                url,
                headers={
                    "User-Agent": "Public-mechanistic-source-collector/1.0",
                    "Accept-Encoding": "identity",
                },
            )
            with urlopen(request, timeout=120) as response:
                raw = response.read(20 * 1024 * 1024 + 1)
                final_url = response.geturl()
                if response.status != 200:
                    raise ValueError("Unexpected HTTP status")
            if len(raw) > 20 * 1024 * 1024:
                raise ValueError("Article response exceeds size limit")
            article, info = parse_article(raw, expected)
            break
        except Exception as error:
            last_error = error
            if attempt == 2:
                raise RuntimeError(
                    f"Could not collect {expected['pmcid']}: {error}"
                ) from error
            time.sleep(2 ** (attempt + 1))
    else:
        raise RuntimeError(str(last_error))

    manifest = {
        "expected_identity": expected,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "requested_url": url,
        "final_url": final_url,
        "article_sha256": sha(raw),
        "article_bytes": len(raw),
        "identity_verified": True,
    }

    with tempfile.TemporaryDirectory(
        prefix=".article_", dir=cache_root
    ) as temporary:
        stage = Path(temporary)
        (stage / "article.xml").write_bytes(raw)
        write_json(stage / "manifest.json", manifest)
        stage.rename(directory)

    return raw, article, info, manifest

def mechanism_claims(info, xml_sha):
    if info["doi"].lower() != "10.1073/pnas.2020599118":
        return []

    definitions = [
        {
            "compound": "ataluren",
            "outcome": "release-factor activity",
            "claim_text": (
                "The authors attribute ataluren-stimulated readthrough "
                "to inhibition of release-factor activity in their "
                "experimental translation work."
            ),
            "quotes": [
                "highly purified in vitro eukaryotic translation system",
                "We find that ataluren stimulation of read-through derives "
                "exclusively from its ability to inhibit release factor activity.",
            ],
        },
        {
            "compound": "G418",
            "outcome": "near-cognate tRNA decoding",
            "claim_text": (
                "The authors report increased functional near-cognate "
                "tRNA mispairing with a premature stop codon for G418, "
                "contrasting its mechanism with ataluren."
            ),
            "quotes": [
                "G418 increases functional near-cognate tRNA mispairing with a PSC",
                "with little if any effect on release factor activity",
            ],
        },
    ]
    rows = []
    for definition in definitions:
        if any(
            norm(quote) not in norm(info["abstract"])
            for quote in definition["quotes"]
        ):
            raise ValueError(
                "Mechanism abstract anchor mismatch: "
                + definition["compound"]
            )

        identity = {
            "source_id": info["source_id"],
            "article_sha256": xml_sha,
            **definition,
        }
        claim_id = "CLM_" + sha(json.dumps(
            identity, sort_keys=True
        ).encode())[:24]

        rows.append({
            "claim_id": claim_id,
            "source_id": info["source_id"],
            "compound": definition["compound"],
            "evidence_class": "purified_in_vitro_translation_experiment",
            "source_gene": "not_disease_gene_specific",
            "population_or_model": (
                "Authors' purified in-vitro eukaryotic translation work; "
                "full-methods adjudication pending"
            ),
            "outcome": definition["outcome"],
            "claim_text": definition["claim_text"],
            "statistics_json": "{}",
            "source_locator": (
                f"sources/{info['pmcid']}.xml:"
                "front/article-meta/abstract"
            ),
            "source_excerpt": " | ".join(definition["quotes"]),
            "limitations": (
                "Mechanistic finding in the studied system. Does not establish "
                "direct binding to a particular protein, eRF1 degradation, "
                "BUB1B rescue, clinical efficacy or current regulatory status."
            ),
            "dependence_group": info["source_id"],
            "extraction_method": "literal_abstract_anchor_with_normalization",
            "review_status": "source_anchors_checked_adjudication_pending",
        })
    return rows

def bib_escape(value):
    mapping = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
        "$": r"\$", "#": r"\#", "_": r"\_",
        "{": r"\{", "}": r"\}",
    }
    return "".join(mapping.get(character, character) for character in value)

def self_test():
    import unittest

    class IdentityTests(unittest.TestCase):
        def expected(self):
            return {
                "pmcid": "PMC1", "pmid": "2",
                "doi": "10.0000/synthetic", "role": "synthetic_test",
            }

        def xml(self, doi="10.0000/synthetic", body=True):
            return (
                '<article><front><article-meta>'
                '<article-id pub-id-type="pmc">1</article-id>'
                '<article-id pub-id-type="pmid">2</article-id>'
                f'<article-id pub-id-type="doi">{doi}</article-id>'
                '<title-group><article-title>Synthetic</article-title>'
                '</title-group></article-meta></front>'
                + ('<body><p>Synthetic body.</p></body>' if body else '')
                + '</article>'
            ).encode()

        def test_identity(self):
            _, info = parse_article(self.xml(), self.expected())
            self.assertEqual(info["source_id"], "PMID2")

        def test_wrong_doi_rejected(self):
            with self.assertRaises(ValueError):
                parse_article(self.xml(doi="10.0000/wrong"), self.expected())

        def test_missing_body_rejected(self):
            with self.assertRaises(ValueError):
                parse_article(self.xml(body=False), self.expected())

        def test_normalization(self):
            self.assertEqual(norm("A  − B"), "A - B")

    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(IdentityTests)
    )
    if not result.wasSuccessful():
        raise SystemExit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--registry")
    parser.add_argument("--cache")
    parser.add_argument("--outdir")
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return
    if not all([args.registry, args.cache, args.outdir]):
        parser.error("--registry, --cache and --outdir are required")

    registry_path = Path(args.registry)
    registry_raw = registry_path.read_bytes()
    registry = json.loads(registry_raw)
    if registry.get("schema_version") != 1:
        raise ValueError("Unsupported registry schema")

    expected_sources = registry["sources"]
    pmcids = [row["pmcid"] for row in expected_sources]
    if len(pmcids) != len(set(pmcids)):
        raise ValueError("Duplicate registry PMCIDs")
    if any(
        not re.fullmatch(r"PMC[0-9]+", row["pmcid"])
        for row in expected_sources
    ):
        raise ValueError("Invalid registry PMCID")

    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    out = Path(args.outdir)
    if out.exists():
        raise ValueError("Output already exists; use a new versioned directory")

    collected = [
        (expected, *cache_article(expected, cache))
        for expected in expected_sources
    ]
    if registry_path.read_bytes() != registry_raw:
        raise ValueError("Registry changed during collection")

    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".mechanisms_", dir=out.parent
    ) as temporary:
        stage = Path(temporary)
        (stage / "sources").mkdir()
        sources, passages, claims, infos, retrievals = [], [], [], [], []
        ris_entries, bib_entries = [], []

        for expected, raw, article, info, retrieval in collected:
            filename = f"sources/{info['pmcid']}.xml"
            (stage / filename).write_bytes(raw)
            xml_sha = sha(raw)

            sources.append({
                "source_id": info["source_id"],
                "source_type": "primary_publication",
                "title": info["title"],
                "authors": info["authors"],
                "year": info["year"],
                "publication": info["publication"],
                "doi": info["doi"],
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{info['pmid']}/",
                "snapshot": filename,
                "snapshot_sha256": xml_sha,
            })
            passages.extend(collect_passages(article, info))
            claims.extend(mechanism_claims(info, xml_sha))
            infos.append(info)
            retrievals.append(retrieval)

            ris = [
                "TY  - JOUR", "ID  - " + info["source_id"],
                "TI  - " + info["title"],
            ]
            ris.extend(
                "AU  - " + name
                for name in info["authors"].split("; ") if name
            )
            for tag, value in [
                ("PY", info["year"]), ("T2", info["publication"]),
                ("VL", info["volume"]), ("IS", info["issue"]),
                ("DO", info["doi"]),
                ("UR", f"https://pubmed.ncbi.nlm.nih.gov/{info['pmid']}/"),
            ]:
                if value:
                    ris.append(f"{tag}  - {value}")
            ris.append("ER  -")
            ris_entries.append("\n".join(ris))

            bib = "@article{" + info["source_id"] + ",\n"
            for field, value in [
                ("title", info["title"]),
                ("author", " and ".join(info["bib_authors"])),
                ("year", info["year"]),
                ("journal", info["publication"]),
                ("volume", info["volume"]),
                ("number", info["issue"]),
                ("doi", info["doi"]),
                ("url", f"https://pubmed.ncbi.nlm.nih.gov/{info['pmid']}/"),
            ]:
                if value:
                    bib += f"  {field} = {{{bib_escape(value)}}},\n"
            bib += "}\n"
            bib_entries.append(bib)

        write_tsv(stage / "sources.tsv", sources, SOURCE_FIELDS)
        write_tsv(stage / "mechanism_claims.tsv", claims, CLAIM_FIELDS)
        write_tsv(
            stage / "review_passages.tsv", passages, PASSAGE_FIELDS
        )
        write_json(stage / "article_metadata.json", infos)
        write_json(stage / "retrieval_provenance.json", retrievals)
        write_json(stage / "source_registry.json", registry)

        write_tsv(
            stage / "review_queue.tsv",
            [{
                "source_id": info["source_id"],
                "role": info["role"],
                "required_review": (
                    "experimental scope and mechanistic inference"
                    if info["role"] == "ataluren_primary_mechanism"
                    else "model inputs, assay normalization, export scale, "
                         "author-context equivalence and functional limitations"
                ),
                "status": "pending",
            } for info in infos],
            ["source_id", "role", "required_review", "status"],
        )

        (stage / "references.ris").write_text(
            "\n\n".join(ris_entries) + "\n", encoding="utf-8"
        )
        (stage / "references.bib").write_text(
            "\n".join(bib_entries), encoding="utf-8"
        )
        (stage / "README.txt").write_text(
            "PRIMARY MECHANISTIC SOURCE AUGMENTATION\n\n"
            "This public-source bundle does not modify an existing case ledger.\n"
            "Sources use PMID citation keys for later deduplication.\n"
            "Mechanism claims are abstract-anchor checked, not adjudicated.\n"
            "Selected full-text passages are review aids, not accepted claims.\n"
            "RTDetective export units are not inferred by this collector.\n"
            "Historical regulatory and broad safety statements are not "
            "converted into current claims.\n"
            "Supplementary files and model/exporter code are not collected.\n"
            "No patient-specific applicability links are created here.\n",
            encoding="utf-8",
        )

        hashes = {
            str(path.relative_to(stage)): sha(path.read_bytes())
            for path in stage.rglob("*") if path.is_file()
        }
        write_json(stage / "manifest.json", {
            "node": "primary_mechanistic_source_collection_v1",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "registry_sha256": sha(registry_raw),
            "implementation_sha256": sha(Path(__file__).read_bytes()),
            "sources": len(sources),
            "mechanism_claims": len(claims),
            "review_passages": len(passages),
            "all_source_identities_verified": True,
            "automatic_scientific_adjudication": False,
            "existing_case_ledger_modified": False,
            "patient_data_read_or_transmitted": False,
            "supplementary_files_or_exporter_code_collected": False,
            "output_sha256": hashes,
        })
        stage.rename(out)

    manifest = json.loads((out / "manifest.json").read_text())
    for filename, expected in manifest["output_sha256"].items():
        if sha((out / filename).read_bytes()) != expected:
            raise ValueError("Published bundle failed checksum verification")

    print("\n===== PRIMARY MECHANISTIC SOURCES COLLECTED =====")
    print("Output:", out)
    print("Sources:", len(sources))
    print("Anchor-checked mechanism claims:", len(claims))
    print("Review passages:", len(passages))
    print("Output checksums: verified")
    for info in infos:
        print(
            info["source_id"], info["pmcid"],
            info["doi"], info["title"], sep="\t"
        )
    print("Existing case ledger remains unchanged.")
    print("RTDetective methods/export-scale interpretation remains pending.")

if __name__ == "__main__":
    main()
