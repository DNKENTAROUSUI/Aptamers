"""
アプタマーデータ読み込みモジュール
==================================
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

KD_CONVERSION = {
    "pm": 1e-3, "nm": 1.0, "um": 1e3, "μm": 1e3, "mm": 1e6, "m": 1e9,
}


def normalize_columns(df):
    df.columns = (
        df.columns.str.strip()
        .str.lower()
        .str.replace(r"\s+", "_", regex=True)
        .str.replace(r"[^\w]", "_", regex=True)
    )
    return df


def filter_protein_targets(df):
    col = "target_type"
    if col not in df.columns:
        warnings.warn(f"\'{col}\' カラムがありません。フィルタリングをスキップします。")
        return df
    before = len(df)
    df = df[df[col].astype(str).str.strip().str.lower() == "protein"].copy()
    logger.info(f"タンパク質ターゲットにフィルタ: {before} → {len(df)} 行")
    return df


def normalize_kd(df):
    if "kd_value" not in df.columns:
        return df
    df["kd_value"] = pd.to_numeric(df["kd_value"], errors="coerce")
    if "kd_unit" in df.columns:
        unit_col = df["kd_unit"].astype(str).str.strip().str.lower()
        factor = unit_col.map(KD_CONVERSION).fillna(1.0)
        df["kd_nm"] = df["kd_value"] * factor
    else:
        df["kd_nm"] = df["kd_value"]
    return df


def deduplicate(df):
    before = len(df)
    df = df.drop_duplicates().reset_index(drop=True)
    logger.info(f"重複除去: {before} → {len(df)} 行")
    return df


def load_and_clean(input_path):
    logger.info(f"入力ファイル読み込み: {input_path}")
    df = safe_read_csv(input_path)
    if df.empty:
        logger.warning("入力データが空です。")
        return df
    df = normalize_columns(df)
    df = filter_protein_targets(df)
    df = normalize_kd(df)
    df = deduplicate(df)

    n_apt = df["aptamer_id"].nunique() if "aptamer_id" in df.columns else len(df)
    n_tgt = df["target_name"].nunique() if "target_name" in df.columns else 0
    n_dna = (df["nucleic_acid_type"].str.upper() == "DNA").sum() if "nucleic_acid_type" in df.columns else 0
    n_rna = (df["nucleic_acid_type"].str.upper() == "RNA").sum() if "nucleic_acid_type" in df.columns else 0
    n_kd_missing = df["kd_nm"].isna().sum() if "kd_nm" in df.columns else 0
    n_source = df["source_database"].value_counts().to_dict() if "source_database" in df.columns else {}

    logger.info("==== データ概要 ====")
    logger.info(f"  アプタマー数     : {n_apt}")
    logger.info(f"  ターゲット数     : {n_tgt}")
    logger.info(f"  DNA              : {n_dna}")
    logger.info(f"  RNA              : {n_rna}")
    logger.info(f"  Kd 欠損          : {n_kd_missing}")
    logger.info(f"  データベース別   : {n_source}")
    logger.info(f"  欠損値合計       : {int(df.isna().sum().sum())}")
    return df


def main():
    parser = argparse.ArgumentParser(description="アプタマーCSV読み込み＆クリーニング")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    df = load_and_clean(args.input)
    safe_write_csv(df, args.output)
    logger.info(f"出力完了: {args.output}  ({len(df)} 行)")


if __name__ == "__main__":
    main()
