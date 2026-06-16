"""
関連解析モジュール
==================
アプタマー構造グループ × タンパク質グループの
クロス集計、オッズ比、Fisher exact / chi-square、多重検定補正、
自然言語サマリー生成を行う。
"""

import argparse
import os
import sys
import warnings
from itertools import product
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from src.analysis.visualization import (
    plot_heatmap, plot_normalized_heatmap, plot_structure_group_counts,
    plot_protein_group_counts, plot_distribution,
)
from src.utils.io import ensure_dir, safe_read_csv, safe_write_csv
from src.utils.logging import setup_logger

logger = setup_logger(__name__)


def prepare_merged_features(aptamer_df, protein_df):
    left_key = "target_name_normalized" if "target_name_normalized" in aptamer_df.columns else "target_name"
    right_key = "target_name_normalized" if "target_name_normalized" in protein_df.columns else "target_name"
    merged = aptamer_df.merge(protein_df, how="left", left_on=left_key, right_on=right_key, suffixes=("_aptamer","_protein"))
    return merged


def compute_crosstab(df, structure_col, protein_group_col):
    return pd.crosstab(df[structure_col], df[protein_group_col])


def odds_ratio_for_pair(df, structure_value, group_value, structure_col, protein_group_col):
    a = ((df[structure_col]==structure_value) & (df[protein_group_col]==group_value)).sum()
    b = ((df[structure_col]==structure_value) & (df[protein_group_col]!=group_value)).sum()
    c = ((df[structure_col]!=structure_value) & (df[protein_group_col]==group_value)).sum()
    d = ((df[structure_col]!=structure_value) & (df[protein_group_col]!=group_value)).sum()
    table = np.array([[a,b],[c,d]])
    table_f = table.astype(float)
    if np.any(table_f==0):
        table_f += 0.5
    or_v = (table_f[0,0]*table_f[1,1])/(table_f[0,1]*table_f[1,0])
    try:
        _, fisher_p = fisher_exact(table)
    except Exception:
        fisher_p = np.nan
    return {
        "structure_group":structure_value, "protein_group":group_value,
        "a":int(a),"b":int(b),"c":int(c),"d":int(d),
        "odds_ratio":float(or_v),
        "fisher_p":float(fisher_p) if pd.notna(fisher_p) else np.nan,
    }


def association_tests(df, structure_col, protein_group_col):
    structures = sorted(df[structure_col].dropna().unique().tolist())
    groups = sorted(df[protein_group_col].dropna().unique().tolist())
    rows = []
    for s, g in product(structures, groups):
        rows.append(odds_ratio_for_pair(df, s, g, structure_col, protein_group_col))
    result = pd.DataFrame(rows)
    if not result.empty and result["fisher_p"].notna().any():
        valid = result["fisher_p"].fillna(1.0).values
        reject, qvals, _, _ = multipletests(valid, method="fdr_bh")
        result["fisher_q"] = qvals
        result["significant_fdr_0_05"] = reject
    else:
        result["fisher_q"] = np.nan
        result["significant_fdr_0_05"] = False
    ctab = pd.crosstab(df[structure_col], df[protein_group_col])
    if ctab.shape[0] >= 2 and ctab.shape[1] >= 2:
        chi2, p, dof, _ = chi2_contingency(ctab)
        logger.info(f"Chi-square: chi2={chi2:.4f}, p={p:.4g}, dof={dof}")
    else:
        logger.warning("chi-square には少なくとも 2x2 が必要")
    return result


def summarize_results(result_df, top_n=10):
    if result_df.empty:
        return "関連解析結果は空でした。"
    sort_df = result_df.sort_values(by=["fisher_q","odds_ratio"], ascending=[True,False], na_position="last")
    lines = ["## 自然言語サマリー","",
             "既知の結合ペア内で、アプタマー構造グループとタンパク質グループの偏りを探索的に調べた。",
             "以下は FDR 補正後 p 値とオッズ比に基づく上位傾向である。",""]
    shown = 0
    for _, row in sort_df.iterrows():
        if shown >= top_n:
            break
        s, g, orv, p, q = row["structure_group"], row["protein_group"], row["odds_ratio"], row["fisher_p"], row["fisher_q"]
        if pd.isna(p):
            continue
        tendency = "多い傾向" if orv > 1 else "少ない傾向"
        lines.append(f"- `{s}` は `{g}` に対して **{tendency}** を示した (odds ratio={orv:.2f}, Fisher p={p:.3g}, FDR q={q:.3g})。")
        shown += 1
    if shown == 0:
        lines.append("- 有効な Fisher 検定結果は十分に得られなかった。")
    lines += ["","※ これは正例データのみを用いた探索的解析であり、因果関係や非結合予測を意味しない。実験的検証が必要。"]
    return "\n".join(lines)


def generate_mvp_report(merged_df, association_df, output_report_path, figures_dir, protein_group_col="charge_group"):
    ensure_dir(os.path.dirname(output_report_path))
    n_pairs = len(merged_df)
    n_targets = merged_df["target_name_normalized"].nunique() if "target_name_normalized" in merged_df.columns else merged_df["target_name"].nunique() if "target_name" in merged_df.columns else 0
    dna_count = (merged_df["nucleic_acid_type"].astype(str).str.upper()=="DNA").sum() if "nucleic_acid_type" in merged_df.columns else 0
    rna_count = (merged_df["nucleic_acid_type"].astype(str).str.upper()=="RNA").sum() if "nucleic_acid_type" in merged_df.columns else 0
    kd_count = merged_df["kd_nm"].notna().sum() if "kd_nm" in merged_df.columns else 0
    structure_counts = merged_df["structure_group"].value_counts().to_dict() if "structure_group" in merged_df.columns else {}
    protein_counts = merged_df[protein_group_col].value_counts().to_dict() if protein_group_col in merged_df.columns else {}
    summary = summarize_results(association_df)
    report = f"""# MVP Report: Aptamer-Protein Structure Association Analysis

## 1. 目的
既知のタンパク質結合アプタマーの配列特徴と簡易構造特徴を抽出し、
タンパク質の物性・グループ分類との関係を探索的に可視化する。

## 2. 使用データ
- アプタマー-タンパク質結合ペア数: {n_pairs}
- ユニークタンパク質数: {n_targets}
- DNA数: {dna_count}
- RNA数: {rna_count}
- Kd値あり: {kd_count}

## 3. データ取得方法
- 入力はローカルCSVを想定
- target_type == protein を対象
- ターゲット名を正規化し、可能なら UniProt 由来のタンパク質配列を利用

## 4. アプタマー特徴量
- 配列長、塩基組成、GC含量
- G4関連特徴 (G-richness, G-skewness, G4Hunter風スコア, motif)
- ステムループ関連特徴 (paired fraction, MFE, stems, hairpins)

## 5. タンパク質特徴量
- アミノ酸配列長、分子量、理論pI、pH 7.4での推定電荷、GRAVY
- 各種アミノ酸群の割合

## 6. 構造グループの分類ルール
- G4_like: G4Hunter風スコア >= 1.0 または G4 motifあり
- stem_loop_like: paired base fraction > 0.3 かつ hairpin loopあり
- G4_and_stem_loop_like: 上記両方
- other_or_uncertain: 上記以外

## 7. タンパク質グループの分類ルール
### 電荷グループ
- positive_charge_like: charge >= 5
- negative_charge_like: charge <= -5
- neutral_charge_like: その間

### 疎水性グループ
- hydrophobic_like: GRAVY >= 0.2
- hydrophilic_like: GRAVY <= -0.5
- intermediate: その間

### pIグループ
- acidic_pI: pI < 6 / neutral_pI: 6-8 / basic_pI: pI > 8

## 8. 主な可視化結果
### 構造グループ件数
{structure_counts}

### タンパク質グループ件数 ({protein_group_col})
{protein_counts}

図:
- `reports/figures/aptamer_structure_counts.png`
- `reports/figures/protein_group_counts.png`
- `reports/figures/structure_vs_protein_group_heatmap.png`

## 9. 統計解析結果
{summary}

## 10. ベースラインモデル結果
MVPでは主に解析パイプラインを優先した。
データ数が十分であれば、後続ステップで RandomForest / LogisticRegression / GradientBoosting を学習可能。

## 11. 限界
- 非結合データを使っていないため、「結合する / しない」の予測ではない
- G4 / stem-loop 判定は簡易推定であり、実験構造を直接表すものではない
- 結果は探索的であり、実験的検証が必要

## 12. 次にやるべきこと
1. より大規模な既知アプタマーデータ統合
2. UniProt マッピング精度向上
3. ViennaRNA / NUPACK ベースの構造推定精密化
4. Negative / decoy データの設計
5. タンパク質グループ予測ベースラインモデルの拡張
"""
    with open(output_report_path, "w", encoding="utf-8") as f:
        f.write(report)
    logger.info(f"レポート保存: {output_report_path}")


def main():
    parser = argparse.ArgumentParser(description="アプタマー構造 × タンパク質グループ関連解析")
    parser.add_argument("--aptamer_features", default="data/processed/aptamer_features.csv")
    parser.add_argument("--protein_features", default="data/processed/protein_features.csv")
    parser.add_argument("--merged_output", default="data/processed/merged_features.csv")
    parser.add_argument("--association_output", default="data/processed/association_results.csv")
    parser.add_argument("--protein_group_col", default="charge_group")
    parser.add_argument("--figures_dir", default="reports/figures")
    parser.add_argument("--report_output", default="reports/mvp_report.md")
    args = parser.parse_args()
    aptamer_df = safe_read_csv(args.aptamer_features)
    protein_df = safe_read_csv(args.protein_features)
    if aptamer_df.empty or protein_df.empty:
        logger.error("入力ファイルが空です。"); return
    merged = prepare_merged_features(aptamer_df, protein_df)
    safe_write_csv(merged, args.merged_output)
    ensure_dir(args.figures_dir)
    if "structure_group" in merged.columns:
        plot_structure_group_counts(merged, os.path.join(args.figures_dir, "aptamer_structure_counts.png"))
    if args.protein_group_col in merged.columns:
        plot_protein_group_counts(merged, args.protein_group_col, os.path.join(args.figures_dir, "protein_group_counts.png"))
    for col in ["sequence_length","gc_content","g_content","mfe","theoretical_pI","estimated_charge_pH7_4","gravy"]:
        if col in merged.columns:
            plot_distribution(merged, col, os.path.join(args.figures_dir, f"{col}_distribution.png"))
    if "structure_group" in merged.columns and args.protein_group_col in merged.columns:
        ctab = compute_crosstab(merged, "structure_group", args.protein_group_col)
        plot_heatmap(ctab, os.path.join(args.figures_dir, "structure_vs_protein_group_heatmap.png"))
        plot_normalized_heatmap(ctab, os.path.join(args.figures_dir, "structure_vs_protein_group_heatmap_normalized.png"))
        assoc = association_tests(merged, "structure_group", args.protein_group_col)
        safe_write_csv(assoc, args.association_output)
        generate_mvp_report(merged, assoc, args.report_output, args.figures_dir, args.protein_group_col)
    else:
        logger.warning("structure_group or protein_group_col not found. Skipping association analysis.")
    logger.info("関連解析完了")

if __name__ == "__main__":
    main()
