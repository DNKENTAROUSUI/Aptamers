"""
UniProt データ取得モジュール
==============================
"""

import argparse
import json
import os
import sys
import time
import warnings
from typing import Any, Dict, List, Optional

import pandas as pd

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from src.utils.logging import setup_logger
from src.utils.io import safe_read_csv, safe_write_csv, ensure_dir

logger = setup_logger(__name__)

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_ORGANISM = "9606"
CACHE_DIR = os.path.join(_PROJECT_ROOT, "data", "raw", ".uniprot_cache")
MAX_RETRIES = 3
RETRY_DELAY = 2.0
REQUEST_INTERVAL = 1.0


def _cache_path(key):
    ensure_dir(CACHE_DIR)
    safe_key = "".join(c if c.isalnum() else "_" for c in key)
    return os.path.join(CACHE_DIR, f"{safe_key}.json")


def _load_cache(key):
    path = _cache_path(key)
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return None


def _save_cache(key, data):
    path = _cache_path(key)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def search_uniprot(protein_name, organism=DEFAULT_ORGANISM):
    cache_key = f"search_{organism}_{protein_name}"
    cached = _load_cache(cache_key)
    if cached is not None:
        logger.info(f"キャッシュヒット: {protein_name}")
        return cached

    try:
        import requests
    except ImportError:
        warnings.warn("requests がインストールされていません。")
        return None

    query = f\'(protein_name:"{protein_name}") AND (organism_id:{organism})\'
    params = {
        "query": query, "format": "json", "size": 1,
        "fields": "accession,gene_names,protein_name,organism_name,sequence,"
                  "keyword,go,cc_subcellular_location",
    }

    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                results = data.get("results", [])
                if not results:
                    logger.warning(f"UniProt 検索結果なし: {protein_name}")
                    empty = {"query": protein_name, "found": False}
                    _save_cache(cache_key, empty)
                    return None
                entry = results[0]
                parsed = _parse_entry(entry, protein_name)
                _save_cache(cache_key, parsed)
                return parsed
            elif resp.status_code == 429:
                time.sleep(RETRY_DELAY * (attempt + 1))
            else:
                logger.warning(f"UniProt API エラー ({resp.status_code}): {protein_name}")
                return None
        except Exception as e:
            logger.warning(f"UniProt API 例外 (attempt {attempt+1}): {e}")
            time.sleep(RETRY_DELAY)
    return None


def _parse_entry(entry, query_name):
    accession = entry.get("primaryAccession", "")
    gene_names_data = entry.get("genes", [])
    gene_name = ""
    if gene_names_data:
        primary = gene_names_data[0].get("geneName", {})
        gene_name = primary.get("value", "")
    protein_desc = entry.get("proteinDescription", {})
    rec_name = protein_desc.get("recommendedName", {})
    prot_name = rec_name.get("fullName", {}).get("value", query_name)
    organism = entry.get("organism", {}).get("scientificName", "")
    seq_data = entry.get("sequence", {})
    sequence = seq_data.get("value", "")
    seq_length = seq_data.get("length", 0)
    mol_weight = seq_data.get("molWeight", 0)
    keywords = [kw.get("name", "") for kw in entry.get("keywords", [])]
    go_terms = []
    xrefs = entry.get("uniProtKBCrossReferences", [])
    for xref in xrefs:
        if xref.get("database") == "GO":
            for p in xref.get("properties", []):
                if p.get("key") == "GoTerm":
                    go_terms.append(p.get("value", ""))
    subcellular = []
    for c in entry.get("comments", []):
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for loc in c.get("subcellularLocations", []):
                val = loc.get("location", {}).get("value", "")
                if val:
                    subcellular.append(val)
    return {
        "query": query_name, "found": True, "accession": accession,
        "gene_name": gene_name, "protein_name": prot_name, "organism": organism,
        "sequence": sequence, "sequence_length": seq_length,
        "molecular_weight": mol_weight, "keywords": keywords,
        "go_terms": go_terms, "subcellular_locations": subcellular,
    }


def batch_fetch_uniprot(target_names, organism=DEFAULT_ORGANISM):
    records = []
    unique_names = list(set(target_names))
    logger.info(f"UniProt 検索開始: {len(unique_names)} 件")
    for i, name in enumerate(unique_names):
        if i > 0:
            time.sleep(REQUEST_INTERVAL)
        logger.info(f"  [{i+1}/{len(unique_names)}] {name}")
        result = search_uniprot(name, organism)
        if result and result.get("found"):
            rec = {
                "target_name": name,
                "uniprot_accession": result.get("accession", ""),
                "gene_name": result.get("gene_name", ""),
                "protein_name_uniprot": result.get("protein_name", ""),
                "organism": result.get("organism", ""),
                "protein_sequence": result.get("sequence", ""),
                "protein_sequence_length": result.get("sequence_length", 0),
                "molecular_weight_da": result.get("molecular_weight", 0),
                "keywords": "|".join(result.get("keywords", [])),
                "go_terms": "|".join(result.get("go_terms", [])),
                "subcellular_locations": "|".join(result.get("subcellular_locations", [])),
            }
        else:
            rec = {
                "target_name": name, "uniprot_accession": "", "gene_name": "",
                "protein_name_uniprot": "", "organism": "", "protein_sequence": "",
                "protein_sequence_length": 0, "molecular_weight_da": 0,
                "keywords": "", "go_terms": "", "subcellular_locations": "",
            }
        records.append(rec)
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="UniProt からタンパク質情報を取得")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--organism", default=DEFAULT_ORGANISM)
    args = parser.parse_args()
    df = safe_read_csv(args.input)
    if "target_name_normalized" in df.columns:
        targets = df["target_name_normalized"].dropna().unique().tolist()
    elif "target_name" in df.columns:
        targets = df["target_name"].dropna().unique().tolist()
    else:
        logger.error("target_name カラムが見つかりません。")
        sys.exit(1)
    result_df = batch_fetch_uniprot(targets, args.organism)
    safe_write_csv(result_df, args.output)
    logger.info(f"出力完了: {args.output} ({len(result_df)} 件)")


if __name__ == "__main__":
    main()
