"""
ターゲットタンパク質名の正規化
================================
"""

import argparse
import os
import sys
import warnings

import pandas as pd

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from src.utils.logging import setup_logger
from src.utils.io import safe_read_csv, safe_write_csv

logger = setup_logger(__name__)

SYNONYM_MAP = {
    "vegf165": "VEGF", "vegf-165": "VEGF", "vegf 165": "VEGF", "vegf": "VEGF",
    "tnf-alpha": "TNF-alpha", "tnf alpha": "TNF-alpha", "tnf-a": "TNF-alpha",
    "tnfa": "TNF-alpha", "tumor necrosis factor alpha": "TNF-alpha",
    "interferon-gamma": "IFN-gamma", "interferon gamma": "IFN-gamma",
    "ifn-gamma": "IFN-gamma", "ifng": "IFN-gamma", "ifn-g": "IFN-gamma",
    "pdgf-bb": "PDGF-BB", "pdgf bb": "PDGF-BB",
    "pdgf-aa": "PDGF-AA", "pdgf aa": "PDGF-AA",
    "ige": "IgE", "immunoglobulin e": "IgE",
    "c-reactive protein": "CRP", "crp": "CRP",
    "egfr": "EGFR", "egf": "EGF", "her2": "HER2", "erbb2": "HER2",
    "psa": "PSA", "prostate specific antigen": "PSA", "prostate-specific antigen": "PSA",
    "human serum albumin": "HSA", "hsa": "HSA", "albumin": "HSA",
    "interleukin-6": "IL-6", "il-6": "IL-6", "il6": "IL-6",
    "ctla-4": "CTLA-4", "ctla4": "CTLA-4",
    "pd-l1": "PD-L1", "pdl1": "PD-L1", "pd-1 ligand": "PD-L1",
    "cd20": "CD20", "fgf2": "FGF2", "fgf-2": "FGF2",
    "basic fibroblast growth factor": "FGF2", "bfgf": "FGF2",
    "lysozyme": "Lysozyme", "thrombin": "Thrombin",
    "streptavidin": "Streptavidin", "insulin": "Insulin",
    "nucleolin": "Nucleolin", "osteopontin": "Osteopontin",
    "muc1": "MUC1", "transferrin": "Transferrin", "fibrinogen": "Fibrinogen",
    "trypsin": "Trypsin", "complement c5": "Complement C5",
    "prion protein": "Prion protein",
    "alpha-fetoprotein": "AFP", "afp": "AFP",
    "angiopoietin-2": "Angiopoietin-2", "ang2": "Angiopoietin-2",
    "cardiac troponin i": "cTnI", "ctni": "cTnI", "troponin i": "cTnI",
    "hemoglobin": "Hemoglobin",
    "neutrophil elastase": "Neutrophil elastase", "elastase": "Neutrophil elastase",
    "pcsk9": "PCSK9",
}


def normalize_target_name(name):
    if not isinstance(name, str) or not name.strip():
        return "Unknown"
    key = name.strip().lower()
    if key in SYNONYM_MAP:
        return SYNONYM_MAP[key]
    return name.strip()


def add_normalized_targets(df):
    if "target_name" not in df.columns:
        warnings.warn("'target_name' カラムがありません。")
        return df
    df = df.copy()
    df["target_name_normalized"] = df["target_name"].apply(normalize_target_name)
    n_changed = (df["target_name"] != df["target_name_normalized"]).sum()
    logger.info(f"ターゲット名正規化: {n_changed} 件が変更されました。")
    return df


def main():
    parser = argparse.ArgumentParser(description="ターゲット名の正規化")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    df = safe_read_csv(args.input)
    df = add_normalized_targets(df)
    safe_write_csv(df, args.output)
    logger.info(f"出力完了: {args.output}")


if __name__ == "__main__":
    main()
