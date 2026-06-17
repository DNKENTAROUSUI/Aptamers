from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import requests

from src.utils.logging import get_logger

LOGGER = get_logger(__name__)
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


@dataclass
class UniProtRecord:
    accession: str | None
    protein_name: str | None
    sequence: str | None
    keywords: str = ""
    go_terms: str = ""
    comments: str = ""


def parse_fasta(path: str | Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current_id: str | None = None
    chunks: list[str] = []
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    records[current_id] = "".join(chunks)
                header = line[1:]
                current_id = header.split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
        if current_id:
            records[current_id] = "".join(chunks)
    return records


def fetch_uniprot_by_target_name(target_name: str, organism: str | None = "9606", sleep_seconds: float = 0.2) -> UniProtRecord:
    query = f'protein_name:"{target_name}"'
    if organism:
        query += f" AND organism_id:{organism}"
    params = {
        "query": query,
        "format": "json",
        "size": 1,
        "fields": "accession,protein_name,sequence,keyword,go,cc_function,cc_subcellular_location",
    }
    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=20)
        response.raise_for_status()
    except requests.RequestException as exc:
        LOGGER.warning("UniProt request failed for %s: %s", target_name, exc)
        return UniProtRecord(None, target_name, None)
    time.sleep(sleep_seconds)
    results = response.json().get("results", [])
    if not results:
        LOGGER.warning("No UniProt match for target: %s", target_name)
        return UniProtRecord(None, target_name, None)
    item = results[0]
    accession = item.get("primaryAccession")
    protein = item.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
    sequence = item.get("sequence", {}).get("value")
    keywords = ";".join(k.get("name", "") for k in item.get("keywords", []))
    go_terms = ";".join(x.get("properties", [{}])[0].get("value", "") for x in item.get("uniProtKBCrossReferences", []) if x.get("database") == "GO")
    comments = ";".join(c.get("commentType", "") + ":" + str(c) for c in item.get("comments", []))
    return UniProtRecord(accession, protein or target_name, sequence, keywords, go_terms, comments)
