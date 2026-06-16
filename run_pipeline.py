"""
簡易パイプライン実行スクリプト
==============================
ローカル CSV 入力から MVP 解析を一通り実行する。
"""

import argparse
import os
import subprocess
import sys


def run(cmd):
    print(f"[RUN] {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")


def main():
    parser = argparse.ArgumentParser(description="Aptamer-protein MVP pipeline runner")
    parser.add_argument("--input", default="data/raw/aptamers_sample.csv")
    parser.add_argument("--skip_uniprot", action="store_true",
                        help="UniProt API 取得をスキップする")
    parser.add_argument("--train_model", action="store_true",
                        help="ベースラインモデルも学習する")
    parser.add_argument("--target", default="charge_group",
                        help="モデルの目的変数 (charge_group / hydrophobicity_group / pI_group)")
    args = parser.parse_args()

    # --- Step 1: データ読み込み・クリーニング ---
    run([
        sys.executable, "-m", "src.data.load_aptamer_data",
        "--input", args.input,
        "--output", "data/processed/aptamer_protein_pairs_clean.csv",
    ])

    # --- Step 2: ターゲット名正規化 ---
    run([
        sys.executable, "-m", "src.data.normalize_targets",
        "--input", "data/processed/aptamer_protein_pairs_clean.csv",
        "--output", "data/processed/aptamer_protein_pairs_clean.csv",
    ])

    # --- Step 3: アプタマー特徴量計算 ---
    run([
        sys.executable, "-m", "src.features.aptamer_features",
        "--input", "data/processed/aptamer_protein_pairs_clean.csv",
        "--output", "data/processed/aptamer_features.csv",
    ])

    # --- Step 4: UniProt 取得 (任意) ---
    if not args.skip_uniprot:
        run([
            sys.executable, "-m", "src.data.fetch_uniprot",
            "--input", "data/processed/aptamer_protein_pairs_clean.csv",
            "--output", "data/processed/protein_features_raw.csv",
        ])
    else:
        print("[INFO] UniProt 取得をスキップしました。")

    # --- Step 5 & 6: タンパク質特徴量 & グループ ---
    raw_prot = "data/processed/protein_features_raw.csv"
    if os.path.exists(raw_prot):
        run([
            sys.executable, "-m", "src.features.protein_features",
            "--input", raw_prot,
            "--output", "data/processed/protein_features_tmp.csv",
        ])
        run([
            sys.executable, "-m", "src.features.protein_grouping",
            "--input", "data/processed/protein_features_tmp.csv",
            "--output", "data/processed/protein_features.csv",
        ])

        # --- Step 7: 関連解析 ---
        run([
            sys.executable, "-m", "src.analysis.association_analysis",
            "--protein_group_col", args.target,
        ])
    else:
        print("[WARN] protein_features_raw.csv が見つかりません。後半パイプラインをスキップ。")

    # --- Step 8: モデル学習 (任意) ---
    if args.train_model and os.path.exists("data/processed/merged_features.csv"):
        run([
            sys.executable, "-m", "src.models.train_baseline",
            "--input", "data/processed/merged_features.csv",
            "--target", args.target,
        ])

    print("\n[DONE] MVP pipeline finished.")


if __name__ == "__main__":
    main()
