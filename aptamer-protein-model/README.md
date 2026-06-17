# Aptamer Protein Model MVP

既知のタンパク質結合アプタマー正例から、アプタマー配列特徴とタンパク質物性・機能グループの関係を探索・可視化するPythonプロジェクトです。

このMVPは「結合する / しない」の厳密な予測ではなく、既知の結合ペアを使った関連分析と、データ数が十分な場合のタンパク質グループ分類ベースラインを対象にします。結果は探索的解析であり、実験的検証が必要です。

## セットアップ

```bash
cd aptamer-protein-model
pip install -r requirements.txt
```

## ローカルCSVで実行する方法

入力CSVは以下のカラムを持つ必要があります。

```csv
aptamer_id,aptamer_name,sequence,nucleic_acid_type,target_name,target_type,kd_value,kd_unit,source_database,reference_doi
```

まず `target_type == protein` の行だけを使います。すぐ試せる拡張seedとして `data/raw/aptamers_literature_seed.csv` を同梱しています。これは55件のアプタマー-タンパク質候補ペアを含み、最初の5件だけの小さい動作確認用CSVより解析向きです。

```bash
python -m src.data.load_aptamer_data \
  --input data/raw/aptamers_literature_seed.csv \
  --output data/processed/aptamer_protein_pairs_clean.csv

python -m src.features.aptamer_features \
  --input data/processed/aptamer_protein_pairs_clean.csv \
  --output data/processed/aptamer_features.csv

python -m src.features.protein_features \
  --input data/processed/aptamer_protein_pairs_clean.csv \
  --output data/processed/protein_features.csv \
  --fetch-uniprot

python -m src.analysis.association_analysis

python -m src.models.train_baseline --target hydrophobicity_group
```

複数CSVを結合して動かす場合:

```bash
python -m src.data.merge_aptamer_sources \
  --inputs data/raw/aptamers_literature_seed.csv data/raw/aptamers_example.csv \
  --output data/processed/aptamer_protein_pairs_clean.csv
```

小さいサンプルデータで動作だけ確認する場合:

```bash
python -m src.data.load_aptamer_data --input data/raw/aptamers_example.csv --output data/processed/aptamer_protein_pairs_clean.csv
python -m src.features.aptamer_features --input data/processed/aptamer_protein_pairs_clean.csv --output data/processed/aptamer_features.csv
python -m src.features.protein_features --input data/processed/aptamer_protein_pairs_clean.csv --fasta data/raw/proteins_example.fasta --output data/processed/protein_features.csv
python -m src.analysis.association_analysis
python -m src.models.train_baseline
```

UniProtを使う場合は、FASTAで見つからない標的だけAPI検索できます。

```bash
python -m src.features.protein_features \
  --input data/processed/aptamer_protein_pairs_clean.csv \
  --output data/processed/protein_features.csv \
  --fetch-uniprot \
  --organism 9606
```

## 出力ファイル

- `data/processed/aptamer_protein_pairs_clean.csv`: protein標的に絞った正規化済みペア
- `data/processed/aptamer_features.csv`: アプタマー配列・G4様・二次構造様特徴
- `data/processed/protein_features.csv`: タンパク質物性・グループ特徴
- `data/processed/merged_features.csv`: 解析・モデル用の結合済み特徴量
- `reports/figures/aptamer_structure_counts.png`: アプタマー構造グループ件数
- `reports/figures/protein_group_counts.png`: タンパク質電荷グループ件数
- `reports/figures/structure_vs_protein_group_heatmap.png`: 構造群 x タンパク質群ヒートマップ
- `reports/mvp_report.md`: 自動生成レポート
- `models/baseline_model.pkl`: データ数が十分な場合に保存されるベースラインモデル。拡張seedでは `hydrophobicity_group` または `pI_group` が学習可能です。

## アプタマー特徴量

基本特徴量として、配列長、A/C/G/T/U割合、GC含量、G含量、C含量、purine/pyrimidine割合、連続G/C最大長、1-mer/2-mer/3-mer頻度、DNA/RNA種別を計算します。

G4様特徴量として、G-richness、G-skewness、G4Hunter風スコア、G-quadruplex motif有無、G-run数、G-run間ループ長を計算します。

二次構造特徴量として、MFE、paired base fraction、stem数、hairpin loop数、最大stem長、最大loop長、dot-bracket構造を計算します。ViennaRNA `RNAfold` が使える場合はそれを使い、ない場合は簡易Nussinov法にフォールバックします。

## 構造グループ分類ルール

閾値は `src/features/structure_grouping.py` にあります。

- `G4_like`: G4 motifに一致、またはG割合、G-run数、G4Hunter風スコアが閾値以上
- `stem_loop_like`: paired base fraction、stem数、hairpin loop数が閾値以上
- `G4_and_stem_loop_like`: 両方を満たす
- `other_or_uncertain`: どちらも満たさない

## タンパク質特徴量と分類

タンパク質配列から、長さ、分子量、理論pI、pH 7.4推定電荷、GRAVY、酸性/塩基性/極性/非極性/芳香族/システイン残基割合を計算します。

閾値とキーワードは `src/features/protein_grouping.py` にあります。

- 電荷: `positive_charge_like`, `negative_charge_like`, `neutral_charge_like`
- 疎水性: `hydrophobic_like`, `hydrophilic_like`, `intermediate`
- pI: `acidic_pI`, `neutral_pI`, `basic_pI`
- 機能・局在: `membrane_protein`, `secreted_or_extracellular`, `cytosolic_or_intracellular`, `enzyme`, `receptor`, `growth_factor_or_cytokine`, `nucleic_acid_binding_protein`, `other_or_unknown`

## MVPの範囲

MVPではローカルCSVから既知の正例を読み込み、アプタマーとタンパク質の特徴量化、クロス集計、ヒートマップ、Fisher exact test/chi-square test、FDR補正、Markdownレポートを作成します。モデル学習はデータ数が十分な場合だけ実行されます。

同梱の `aptamers_literature_seed.csv` は解析パイプライン検証用のseedデータです。引用・単位・配列修飾の表記は今後の外部DB連携で精査する前提で、最終的な研究利用では元論文または公式DBで再確認してください。

## 今後の拡張案

- AptaDB、UTexas Aptamer Database、Ribocentre、Aptamer Base、Apta-Indexの取得コネクタ追加
- 標的名からUniProt accessionへの正規化強化
- 実験条件、Kd単位、配列切り出しの標準化
- RNAfold/NUPACKによる二次構造特徴の本格化
- 負例または候補タンパク質集合を設計した予測タスク化
- 未知アプタマー入力CLIと説明可能性レポートの追加
