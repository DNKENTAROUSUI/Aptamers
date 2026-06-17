from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests

from src.analysis.visualization import save_bar_counts, save_heatmap, save_hist
from src.utils.io import read_csv, write_csv
from src.utils.logging import get_logger

LOGGER = get_logger(__name__)
PROTEIN_GROUP_COLUMNS = ["charge_group", "hydrophobicity_group", "pI_group", "functional_group"]


def merge_features(pairs: pd.DataFrame, aptamer: pd.DataFrame, protein: pd.DataFrame) -> pd.DataFrame:
    apt_cols = [c for c in aptamer.columns if c not in {"target_name"}]
    merge_key = ["aptamer_id"] if "aptamer_id" in pairs.columns and "aptamer_id" in aptamer.columns else ["sequence"]
    merged = pairs.merge(aptamer[apt_cols], on=merge_key, how="left", suffixes=("", "_aptamer_feature"))
    merged = merged.merge(protein, on="target_name", how="left")
    return merged


def data_overview(pairs: pd.DataFrame) -> dict[str, object]:
    return {
        "aptamer_count": int(pairs["aptamer_id"].nunique()),
        "target_count": int(pairs["target_name"].nunique()),
        "nucleic_acid_type_counts": pairs["nucleic_acid_type"].fillna("unknown").value_counts().to_dict(),
        "kd_available_count": int(pd.to_numeric(pairs["kd_value"], errors="coerce").notna().sum()),
        "source_database_counts": pairs["source_database"].fillna("unknown").value_counts().to_dict(),
        "missing_values": pairs.isna().sum().to_dict(),
    }


def association_tests(df: pd.DataFrame, protein_group_col: str) -> pd.DataFrame:
    rows = []
    subset = df[["aptamer_structure_group", protein_group_col]].dropna()
    subset = subset[(subset["aptamer_structure_group"] != "") & (subset[protein_group_col] != "unknown")]
    if subset.empty:
        return pd.DataFrame()
    for apt_group in sorted(subset["aptamer_structure_group"].unique()):
        for prot_group in sorted(subset[protein_group_col].unique()):
            a = int(((subset["aptamer_structure_group"] == apt_group) & (subset[protein_group_col] == prot_group)).sum())
            b = int(((subset["aptamer_structure_group"] == apt_group) & (subset[protein_group_col] != prot_group)).sum())
            c = int(((subset["aptamer_structure_group"] != apt_group) & (subset[protein_group_col] == prot_group)).sum())
            d = int(((subset["aptamer_structure_group"] != apt_group) & (subset[protein_group_col] != prot_group)).sum())
            table = np.array([[a, b], [c, d]])
            odds_ratio, p_value = fisher_exact(table) if table.sum() else (np.nan, np.nan)
            rows.append(
                {
                    "protein_group_column": protein_group_col,
                    "aptamer_structure_group": apt_group,
                    "protein_group": prot_group,
                    "count_in_cell": a,
                    "odds_ratio": odds_ratio,
                    "p_value": p_value,
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["p_value_fdr_bh"] = multipletests(result["p_value"].fillna(1.0), method="fdr_bh")[1]
    return result


def chi_square_summary(df: pd.DataFrame, protein_group_col: str) -> dict[str, object]:
    table = pd.crosstab(df["aptamer_structure_group"], df[protein_group_col])
    if table.shape[0] < 2 or table.shape[1] < 2:
        return {"protein_group_column": protein_group_col, "chi2": np.nan, "p_value": np.nan, "dof": np.nan}
    chi2, p, dof, _ = chi2_contingency(table)
    return {"protein_group_column": protein_group_col, "chi2": chi2, "p_value": p, "dof": dof}


def natural_language_summary(stats: pd.DataFrame) -> list[str]:
    if stats.empty:
        return ["統計検定に十分なカテゴリデータがありませんでした。"]
    significant = stats[(stats["p_value_fdr_bh"] < 0.05) & (stats["odds_ratio"] > 1)]
    if significant.empty:
        return ["FDR補正後p<0.05で明確に過剰出現する構造群とタンパク質群の組み合わせは見つかりませんでした。"]
    lines = []
    for _, row in significant.sort_values("p_value_fdr_bh").head(10).iterrows():
        lines.append(
            f"{row['aptamer_structure_group']} は {row['protein_group']} ({row['protein_group_column']}) に多い傾向があります "
            f"(odds ratio={row['odds_ratio']:.2f}, FDR p={row['p_value_fdr_bh']:.3g})。"
        )
    return lines


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return "該当なし"
    headers = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for _, row in df.iterrows():
        values = []
        for col in df.columns:
            value = row[col]
            if isinstance(value, float):
                values.append(f"{value:.4g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_report(
    output: str | Path,
    overview: dict[str, object],
    chi2: pd.DataFrame,
    fisher: pd.DataFrame,
    model_summary: str = "モデル学習は未実行、またはデータ数不足のためスキップされました。",
) -> None:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    summary_lines = natural_language_summary(fisher)
    text = f"""# Aptamer-Protein MVP Report

## 1. 目的
既知のタンパク質結合アプタマー正例から、アプタマー配列・構造様特徴とタンパク質物性・機能グループの関係を探索的に解析する。

## 2. 使用データ
- アプタマー数: {overview.get('aptamer_count')}
- タンパク質標的数: {overview.get('target_count')}
- DNA/RNA比率: {overview.get('nucleic_acid_type_counts')}
- Kd値あり: {overview.get('kd_available_count')}
- データベース別件数: {overview.get('source_database_counts')}

## 3. データ取得方法
MVPではローカルCSVを優先する。タンパク質配列はローカルFASTA、または `--fetch-uniprot` 指定時にUniProt REST APIから取得する。

## 4. アプタマー特徴量
配列長、塩基割合、GC含量、purine/pyrimidine割合、G/C連続長、1-3 mer頻度、G4様スコア、G-run、MFE、paired base fraction、stem/loop特徴を計算した。

## 5. タンパク質特徴量
配列長、分子量、理論pI、pH 7.4推定電荷、GRAVY、残基カテゴリ割合を計算した。

## 6. 構造グループの分類ルール
`G4_like` はG4モチーフ、またはG割合・G-run数・G4Hunter風スコアの閾値で判定する。`stem_loop_like` はpaired base fraction、stem数、hairpin loop数で判定する。両方を満たす場合は `G4_and_stem_loop_like`、どちらでもない場合は `other_or_uncertain` とする。

## 7. タンパク質グループの分類ルール
電荷はpH 7.4推定電荷、疎水性はGRAVY、pIは理論pIで分類する。機能・局在はUniProt annotation/keyword/GO由来テキストのキーワードマッチで分類する。

## 8. 主な可視化結果
図は `reports/figures/` に保存した。主要出力はアプタマー構造群件数、タンパク質グループ件数、構造群とタンパク質グループのヒートマップである。

## 9. 統計解析結果
Chi-square summary:

{dataframe_to_markdown(chi2)}

Fisher exact test summary:

{dataframe_to_markdown(fisher.sort_values('p_value_fdr_bh').head(20) if not fisher.empty else fisher)}

自然言語サマリー:
{chr(10).join('- ' + line for line in summary_lines)}

## 10. ベースラインモデル結果
{model_summary}

## 11. 限界
この結果は既知の正例データに基づく探索的解析であり、非結合を直接予測するものではない。外部DB由来データには選択バイアス、表記揺れ、Kd条件差、配列切り出し差がある。二次構造は予測であり、実験的検証が必要である。

## 12. 次にやるべきこと
データソースを拡充し、標的名の正規化とUniProt accession対応を改善する。測定条件を標準化し、負例または候補集合を設計したうえで、独立テストセットによる予測性能評価を行う。
"""
    Path(output).write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", default="data/processed/aptamer_protein_pairs_clean.csv")
    parser.add_argument("--aptamer-features", default="data/processed/aptamer_features.csv")
    parser.add_argument("--protein-features", default="data/processed/protein_features.csv")
    parser.add_argument("--merged-output", default="data/processed/merged_features.csv")
    parser.add_argument("--stats-output", default="data/processed/association_stats.csv")
    parser.add_argument("--report-output", default="reports/mvp_report.md")
    parser.add_argument("--figures-dir", default="reports/figures")
    parser.add_argument("--model-summary", default="data/processed/baseline_model_summary.md")
    args = parser.parse_args()

    pairs = read_csv(args.pairs)
    aptamer = read_csv(args.aptamer_features)
    protein = read_csv(args.protein_features)
    merged = merge_features(pairs, aptamer, protein)
    write_csv(merged, args.merged_output)

    fig_dir = Path(args.figures_dir)
    save_bar_counts(merged["aptamer_structure_group"], "Aptamer structure groups", fig_dir / "aptamer_structure_counts.png")
    save_bar_counts(merged["charge_group"], "Protein charge groups", fig_dir / "protein_group_counts.png")
    save_heatmap(pd.crosstab(merged["aptamer_structure_group"], merged["charge_group"]), "Aptamer structure vs protein charge group", fig_dir / "structure_vs_protein_group_heatmap.png")
    save_hist(merged["sequence_length"], "Aptamer sequence length distribution", fig_dir / "aptamer_sequence_length_distribution.png")
    save_hist(merged["GC_content"], "Aptamer GC content distribution", fig_dir / "aptamer_gc_content_distribution.png")
    save_hist(merged["G_content"], "Aptamer G content distribution", fig_dir / "aptamer_g_content_distribution.png")
    save_hist(merged["MFE"], "Predicted MFE distribution", fig_dir / "aptamer_mfe_distribution.png")
    save_hist(merged["theoretical_pI"], "Protein pI distribution", fig_dir / "protein_pi_distribution.png")
    save_hist(merged["charge_at_pH_7_4"], "Protein charge distribution", fig_dir / "protein_charge_distribution.png")
    save_hist(merged["GRAVY"], "Protein GRAVY distribution", fig_dir / "protein_gravy_distribution.png")

    fisher_tables = [association_tests(merged, col) for col in PROTEIN_GROUP_COLUMNS if col in merged.columns]
    fisher = pd.concat([t for t in fisher_tables if not t.empty], ignore_index=True) if fisher_tables else pd.DataFrame()
    chi2 = pd.DataFrame([chi_square_summary(merged, col) for col in PROTEIN_GROUP_COLUMNS if col in merged.columns])
    write_csv(fisher, args.stats_output)
    model_summary = "モデル学習は未実行、またはデータ数不足のためスキップされました。"
    if Path(args.model_summary).exists():
        model_summary = Path(args.model_summary).read_text(encoding="utf-8")
    write_report(args.report_output, data_overview(pairs), chi2, fisher, model_summary=model_summary)
    LOGGER.info("Wrote merged data, figures, association stats, and report.")


if __name__ == "__main__":
    main()
