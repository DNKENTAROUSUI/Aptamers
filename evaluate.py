"""
タンパク質物性特徴量計算モジュール
================================
"""

import argparse
import math
import os
import sys
import warnings
from typing import Dict

import numpy as np
import pandas as pd

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from src.utils.logging import setup_logger
from src.utils.io import safe_read_csv, safe_write_csv

logger = setup_logger(__name__)

PH_FOR_CHARGE = 7.4

ACIDIC_RESIDUES  = set("DE")
BASIC_RESIDUES   = set("RKH")
POLAR_RESIDUES   = set("STNQYC")
NONPOLAR_RESIDUES= set("AVILMFWPG")
AROMATIC_RESIDUES= set("FWY")

HYDROPATHY = {
    "A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,
    "G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,
    "P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2,
}

AA_MW = {
    "A":89.09,"R":174.20,"N":132.12,"D":133.10,"C":121.15,"Q":146.15,
    "E":147.13,"G":75.07,"H":155.16,"I":131.17,"L":131.17,"K":146.19,
    "M":149.21,"F":165.19,"P":115.13,"S":105.09,"T":119.12,"W":204.23,
    "Y":181.19,"V":117.15,
}

PKA_NTERM = 9.69
PKA_CTERM = 2.34
PKA_SIDE = {"C":8.33,"D":3.86,"E":4.25,"H":6.00,"K":10.53,"R":12.48,"Y":10.07}


def clean_protein_sequence(sequence: str) -> str:
    if not isinstance(sequence, str):
        return ""
    seq = sequence.strip().upper()
    valid = set(AA_MW.keys())
    return "".join([aa for aa in seq if aa in valid])


def aa_fraction(seq: str, residue_set: set) -> float:
    if not seq:
        return 0.0
    return sum(1 for aa in seq if aa in residue_set) / len(seq)


def manual_molecular_weight(seq: str) -> float:
    if not seq:
        return 0.0
    mw = sum(AA_MW.get(aa, 0.0) for aa in seq)
    water_loss = 18.015 * (len(seq) - 1) if len(seq) > 1 else 0.0
    return mw - water_loss


def manual_gravy(seq: str) -> float:
    if not seq:
        return 0.0
    vals = [HYDROPATHY.get(aa, 0.0) for aa in seq]
    return float(np.mean(vals)) if vals else 0.0


def charge_of_group(pH: float, pKa: float, positive: bool) -> float:
    if positive:
        return 1.0 / (1.0 + 10 ** (pH - pKa))
    else:
        return -1.0 / (1.0 + 10 ** (pKa - pH))


def estimate_net_charge(seq: str, pH: float = PH_FOR_CHARGE) -> float:
    if not seq:
        return 0.0
    charge = 0.0
    charge += charge_of_group(pH, PKA_NTERM, positive=True)
    charge += charge_of_group(pH, PKA_CTERM, positive=False)
    counts = {aa: seq.count(aa) for aa in set(seq)}
    for aa in ("K","R","H"):
        charge += counts.get(aa, 0) * charge_of_group(pH, PKA_SIDE[aa], positive=True)
    for aa in ("D","E","C","Y"):
        charge += counts.get(aa, 0) * charge_of_group(pH, PKA_SIDE[aa], positive=False)
    return charge


def estimate_pI(seq: str) -> float:
    if not seq:
        return 0.0
    low, high = 0.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2
        c = estimate_net_charge(seq, pH=mid)
        if c > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2, 3)


def compute_protein_features(sequence: str) -> Dict:
    seq = clean_protein_sequence(sequence)
    empty = {
        "protein_sequence_length":0,"molecular_weight_da":0.0,"theoretical_pI":0.0,
        "estimated_charge_pH7_4":0.0,"gravy":0.0,"acidic_fraction":0.0,
        "basic_fraction":0.0,"polar_fraction":0.0,"nonpolar_fraction":0.0,
        "aromatic_fraction":0.0,"cysteine_fraction":0.0,
    }
    if not seq:
        return empty
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        analysis = ProteinAnalysis(seq)
        mw = analysis.molecular_weight()
        pI = analysis.isoelectric_point()
        gravy = analysis.gravy()
    except Exception:
        mw = manual_molecular_weight(seq)
        pI = estimate_pI(seq)
        gravy = manual_gravy(seq)
    charge_7_4 = estimate_net_charge(seq, PH_FOR_CHARGE)
    return {
        "protein_sequence_length": len(seq),
        "molecular_weight_da": round(float(mw), 3),
        "theoretical_pI": round(float(pI), 3),
        "estimated_charge_pH7_4": round(float(charge_7_4), 3),
        "gravy": round(float(gravy), 4),
        "acidic_fraction": round(aa_fraction(seq, ACIDIC_RESIDUES), 4),
        "basic_fraction": round(aa_fraction(seq, BASIC_RESIDUES), 4),
        "polar_fraction": round(aa_fraction(seq, POLAR_RESIDUES), 4),
        "nonpolar_fraction": round(aa_fraction(seq, NONPOLAR_RESIDUES), 4),
        "aromatic_fraction": round(aa_fraction(seq, AROMATIC_RESIDUES), 4),
        "cysteine_fraction": round(seq.count("C") / len(seq), 4),
    }


def compute_all_protein_features(df: pd.DataFrame) -> pd.DataFrame:
    seq_col = None
    for c in ["protein_sequence","sequence","aa_sequence"]:
        if c in df.columns:
            seq_col = c
            break
    if seq_col is None:
        warnings.warn("protein sequence column が見つかりません。")
        return df
    logger.info(f"タンパク質特徴量計算開始: {len(df)} 行")
    feats = []
    for idx, row in df.iterrows():
        seq = row.get(seq_col, "")
        try:
            feats.append(compute_protein_features(seq))
        except Exception as e:
            warnings.warn(f"行 {idx} エラー: {e}")
            feats.append(compute_protein_features(""))
    feat_df = pd.DataFrame(feats)
    keep_cols = [c for c in [
        "target_name","target_name_normalized","uniprot_accession",
        "gene_name","protein_name_uniprot","organism",
        seq_col,"keywords","go_terms","subcellular_locations"
    ] if c in df.columns]
    return pd.concat([df[keep_cols].reset_index(drop=True), feat_df], axis=1)


def main():
    parser = argparse.ArgumentParser(description="タンパク質特徴量計算")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    df = safe_read_csv(args.input)
    result = compute_all_protein_features(df)
    safe_write_csv(result, args.output)
    logger.info(f"出力完了: {args.output} ({len(result)} 行)")

if __name__ == "__main__":
    main()
