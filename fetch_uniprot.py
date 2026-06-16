"""
可視化モジュール
================
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from src.utils.io import ensure_dir
from src.utils.logging import setup_logger

logger = setup_logger(__name__)
sns.set_context("talk")


def plot_structure_group_counts(df, output_path):
    ensure_dir(os.path.dirname(output_path))
    counts = df["structure_group"].value_counts()
    plt.figure(figsize=(10,6))
    sns.barplot(x=counts.index, y=counts.values, palette="viridis")
    plt.xticks(rotation=20, ha="right")
    plt.title("Aptamer structure group counts")
    plt.ylabel("Count"); plt.xlabel("Structure group")
    plt.tight_layout(); plt.savefig(output_path, dpi=200); plt.close()
    logger.info(f"保存: {output_path}")


def plot_protein_group_counts(df, group_col, output_path):
    ensure_dir(os.path.dirname(output_path))
    counts = df[group_col].value_counts()
    plt.figure(figsize=(10,6))
    sns.barplot(x=counts.index, y=counts.values, palette="magma")
    plt.xticks(rotation=25, ha="right")
    plt.title(f"Protein {group_col} counts")
    plt.ylabel("Count"); plt.xlabel(group_col)
    plt.tight_layout(); plt.savefig(output_path, dpi=200); plt.close()
    logger.info(f"保存: {output_path}")


def plot_distribution(df, col, output_path, bins=20):
    ensure_dir(os.path.dirname(output_path))
    plt.figure(figsize=(8,6))
    sns.histplot(df[col].dropna(), bins=bins, kde=True)
    plt.title(f"{col} distribution"); plt.xlabel(col)
    plt.tight_layout(); plt.savefig(output_path, dpi=200); plt.close()
    logger.info(f"保存: {output_path}")


def plot_heatmap(crosstab_df, output_path, annot=True):
    ensure_dir(os.path.dirname(output_path))
    plt.figure(figsize=(10,7))
    sns.heatmap(crosstab_df, annot=annot, fmt=".0f", cmap="YlGnBu")
    plt.title("Aptamer structure vs protein group")
    plt.xlabel("Protein group"); plt.ylabel("Aptamer structure group")
    plt.tight_layout(); plt.savefig(output_path, dpi=220); plt.close()
    logger.info(f"保存: {output_path}")


def plot_normalized_heatmap(crosstab_df, output_path, annot=True):
    ensure_dir(os.path.dirname(output_path))
    row_norm = crosstab_df.div(crosstab_df.sum(axis=1).replace(0,1), axis=0)
    plt.figure(figsize=(10,7))
    sns.heatmap(row_norm, annot=annot, fmt=".2f", cmap="coolwarm", center=0.5)
    plt.title("Normalized structure vs protein group heatmap")
    plt.xlabel("Protein group"); plt.ylabel("Aptamer structure group")
    plt.tight_layout(); plt.savefig(output_path, dpi=220); plt.close()
    logger.info(f"保存: {output_path}")
