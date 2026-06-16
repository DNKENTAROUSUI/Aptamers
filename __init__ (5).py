"""
アプタマー構造グループ分類モジュール
======================================
"""

import os
import sys
import warnings

import pandas as pd

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from src.utils.logging import setup_logger

logger = setup_logger(__name__)

GROUP_G4 = "G4_like"
GROUP_STEM_LOOP = "stem_loop_like"
GROUP_BOTH = "G4_and_stem_loop_like"
GROUP_OTHER = "other_or_uncertain"
ALL_STRUCTURE_GROUPS = [GROUP_G4, GROUP_STEM_LOOP, GROUP_BOTH, GROUP_OTHER]


def classify_structure_group(is_g4_like, is_stem_loop_like):
    if is_g4_like and is_stem_loop_like:
        return GROUP_BOTH
    elif is_g4_like:
        return GROUP_G4
    elif is_stem_loop_like:
        return GROUP_STEM_LOOP
    return GROUP_OTHER


def assign_structure_groups(df):
    required = {"is_g4_like", "is_stem_loop_like"}
    if not required.issubset(df.columns):
        warnings.warn(f"必要なカラムがありません: {required - set(df.columns)}")
        return df
    df = df.copy()
    df["structure_group"] = df.apply(
        lambda r: classify_structure_group(bool(r["is_g4_like"]), bool(r["is_stem_loop_like"])),
        axis=1,
    )
    return df


def structure_group_summary(df):
    if "structure_group" not in df.columns:
        warnings.warn("'structure_group' カラムがありません。")
        return pd.Series(dtype=int)
    return df["structure_group"].value_counts().reindex(ALL_STRUCTURE_GROUPS, fill_value=0)
