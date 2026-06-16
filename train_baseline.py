"""
タンパク質グループ分類モジュール
==============================
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

CHARGE_POS_THRESHOLD = 5.0
CHARGE_NEG_THRESHOLD = -5.0
GRAVY_HYDROPHOBIC_THRESHOLD = 0.2
GRAVY_HYDROPHILIC_THRESHOLD = -0.5
PI_ACIDIC_THRESHOLD = 6.0
PI_BASIC_THRESHOLD = 8.0

MEMBRANE_KEYWORDS = {"membrane","transmembrane","cell membrane","plasma membrane"}
SECRETED_KEYWORDS = {"secreted","extracellular","extracellular region","extracellular space"}
INTRACELLULAR_KEYWORDS = {"cytoplasm","cytosol","nucleus","intracellular","nucleoplasm"}
ENZYME_KEYWORDS = {"enzyme","catalytic activity","hydrolase","kinase","protease","transferase","oxidoreductase","ligase","isomerase","lyase"}
RECEPTOR_KEYWORDS = {"receptor","signal transducer","transmembrane receptor"}
CYTOKINE_KEYWORDS = {"cytokine","growth factor","interleukin","interferon","tumor necrosis factor","chemokine"}
NUCLEIC_ACID_BINDING_KEYWORDS = {"dna-binding","rna-binding","nucleic acid binding","nucleolin"}


def classify_charge_group(charge):
    if pd.isna(charge): return "neutral_charge_like"
    if charge >= CHARGE_POS_THRESHOLD: return "positive_charge_like"
    if charge <= CHARGE_NEG_THRESHOLD: return "negative_charge_like"
    return "neutral_charge_like"

def classify_hydrophobicity_group(gravy):
    if pd.isna(gravy): return "intermediate"
    if gravy >= GRAVY_HYDROPHOBIC_THRESHOLD: return "hydrophobic_like"
    if gravy <= GRAVY_HYDROPHILIC_THRESHOLD: return "hydrophilic_like"
    return "intermediate"

def classify_pI_group(pI):
    if pd.isna(pI): return "neutral_pI"
    if pI < PI_ACIDIC_THRESHOLD: return "acidic_pI"
    if pI > PI_BASIC_THRESHOLD: return "basic_pI"
    return "neutral_pI"

def contains_any(text, keywords):
    text = (text or "").lower()
    return any(k in text for k in keywords)

def classify_functional_group(row):
    text = " ".join([
        str(row.get("protein_name_uniprot","")),str(row.get("keywords","")),
        str(row.get("go_terms","")),str(row.get("subcellular_locations","")),
        str(row.get("target_name_normalized","")),str(row.get("target_name","")),
    ]).lower()
    if contains_any(text, MEMBRANE_KEYWORDS): return "membrane_protein"
    if contains_any(text, SECRETED_KEYWORDS): return "secreted_or_extracellular"
    if contains_any(text, RECEPTOR_KEYWORDS): return "receptor"
    if contains_any(text, CYTOKINE_KEYWORDS): return "growth_factor_or_cytokine"
    if contains_any(text, ENZYME_KEYWORDS): return "enzyme"
    if contains_any(text, NUCLEIC_ACID_BINDING_KEYWORDS): return "nucleic_acid_binding_protein"
    if contains_any(text, INTRACELLULAR_KEYWORDS): return "cytosolic_or_intracellular"
    return "other_or_unknown"


def add_protein_groups(df):
    df = df.copy()
    df["charge_group"] = df.get("estimated_charge_pH7_4", pd.Series([0]*len(df))).apply(classify_charge_group)
    df["hydrophobicity_group"] = df.get("gravy", pd.Series([0]*len(df))).apply(classify_hydrophobicity_group)
    df["pI_group"] = df.get("theoretical_pI", pd.Series([7]*len(df))).apply(classify_pI_group)
    df["functional_group"] = df.apply(classify_functional_group, axis=1)
    return df


def main():
    parser = argparse.ArgumentParser(description="タンパク質グループ分類")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    df = safe_read_csv(args.input)
    result = add_protein_groups(df)
    safe_write_csv(result, args.output)
    logger.info(f"出力完了: {args.output}")

if __name__ == "__main__":
    main()
