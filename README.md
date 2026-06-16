# aptamer-protein-model

アプタマー配列特徴とタンパク質物性の関係を探索的に解析し、
**「どのような構造特徴を持つアプタマーが、どのような性質のタンパク質に偏りやすいか」**
を可視化するための Python プロジェクトです。

> **注意**: この結果は探索的解析であり、実験的検証が必要です。

---

## 1. このプロジェクトの目的

MVP の目的は、**未知タンパク質への厳密な結合予測**ではなく、
**既知のアプタマー–タンパク質結合データから、構造とタンパク質性質の関係を見える化すること**です。

1. 既知のタンパク質結合アプタマーをローカル CSV から読み込む
2. アプタマー配列特徴量を計算する（GC含量、k-mer、G4スコア、二次構造推定…）
3. G4_like / stem_loop_like を推定・分類する
4. タンパク質配列から物性特徴量（pI, 電荷, GRAVY, アミノ酸組成）を計算する
5. タンパク質を電荷・疎水性・pI・機能で分類する
6. アプタマー構造グループ × タンパク質グループをクロス集計しヒートマップを出す
7. Fisher exact test / chi-square / odds ratio で偏りを評価する
8. 簡易レポート `reports/mvp_report.md` を自動生成する

---

## 2. MVP の範囲

| 項目 | 内容 |
|------|------|
| 入力 | アプタマー配列、標的タンパク質名、Kd、DOI などを含む CSV |
| データ対象 | タンパク質に結合すると報告されたアプタマーのみ（正例のみ） |
| アプタマー分類 | `G4_like` / `stem_loop_like` / `G4_and_stem_loop_like` / `other_or_uncertain` |
| タンパク質分類 | 電荷、pI、疎水性、機能分類 |
| 出力 | クロス集計、ヒートマップ、統計検定、簡易レポート |
| モデル | 余裕があれば RandomForest / LogisticRegression / GradientBoosting のベースライン |

---

## 3. データの取得方法

### 方法A: ローカル CSV を使う（推奨）

以下のカラムを含む CSV を用意します。

```csv
aptamer_id,aptamer_name,sequence,nucleic_acid_type,target_name,target_type,kd_value,kd_unit,source_database,reference_doi
```

- `target_type == protein` のみが解析対象
- `kd_unit` は pM / nM / uM / mM に対応（自動で nM に変換）
- Kd 欠損があっても動作

サンプルデータ: `data/raw/aptamers_sample.csv`（55 エントリ）

### 方法B: UniProt API を使う

```bash
python -m src.data.fetch_uniprot \
  --input data/processed/aptamer_protein_pairs_clean.csv \
  --output data/processed/protein_features_raw.csv
```

ただし API は不安定なため、**まずローカル CSV で動くことを優先**しています。

---

## 4. 環境構築方法

```bash
git clone <your-repo-url>
cd aptamer-protein-model

python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate

pip install -r requirements.txt
```

### 依存パッケージ

- pandas, numpy, biopython, scikit-learn
- matplotlib, seaborn, scipy, statsmodels
- requests, tqdm, joblib

---

## 5. 実行コマンド

### 一括実行（推奨）

```bash
# UniProt API を使う場合
python run_pipeline.py --input data/raw/aptamers_sample.csv

# UniProt API をスキップする場合
python run_pipeline.py --input data/raw/aptamers_sample.csv --skip_uniprot

# モデル学習も含める場合
python run_pipeline.py --input data/raw/aptamers_sample.csv --train_model --target charge_group
```

### ステップバイステップ実行

```bash
# Step 1: データ読み込み・クリーニング
python -m src.data.load_aptamer_data \
  --input data/raw/aptamers_sample.csv \
  --output data/processed/aptamer_protein_pairs_clean.csv

# Step 2: ターゲット名正規化
python -m src.data.normalize_targets \
  --input data/processed/aptamer_protein_pairs_clean.csv \
  --output data/processed/aptamer_protein_pairs_clean.csv

# Step 3: アプタマー特徴量計算
python -m src.features.aptamer_features \
  --input data/processed/aptamer_protein_pairs_clean.csv \
  --output data/processed/aptamer_features.csv

# Step 4: UniProt 取得（任意）
python -m src.data.fetch_uniprot \
  --input data/processed/aptamer_protein_pairs_clean.csv \
  --output data/processed/protein_features_raw.csv

# Step 5: タンパク質特徴量計算
python -m src.features.protein_features \
  --input data/processed/protein_features_raw.csv \
  --output data/processed/protein_features_tmp.csv

# Step 6: タンパク質グループ分類
python -m src.features.protein_grouping \
  --input data/processed/protein_features_tmp.csv \
  --output data/processed/protein_features.csv

# Step 7: 関連解析・可視化・レポート生成
python -m src.analysis.association_analysis

# Step 8: モデル学習（任意）
python -m src.models.train_baseline \
  --input data/processed/merged_features.csv \
  --target charge_group
```

---

## 6. テストの実行

```bash
pytest tests/ -v
```

---

## 7. 出力ファイルの説明

### データ

| ファイル | 内容 |
|----------|------|
| `data/processed/aptamer_protein_pairs_clean.csv` | クリーニング済みペアデータ |
| `data/processed/aptamer_features.csv` | アプタマー特徴量 |
| `data/processed/protein_features.csv` | タンパク質特徴量＋グループ |
| `data/processed/merged_features.csv` | 統合特徴量 |
| `data/processed/association_results.csv` | 統計検定結果 |

### 図

| ファイル | 内容 |
|----------|------|
| `reports/figures/aptamer_structure_counts.png` | 構造グループ件数 |
| `reports/figures/protein_group_counts.png` | タンパク質グループ件数 |
| `reports/figures/structure_vs_protein_group_heatmap.png` | クロス集計ヒートマップ |
| `reports/figures/structure_vs_protein_group_heatmap_normalized.png` | 正規化ヒートマップ |
| `reports/figures/*_distribution.png` | 各種分布図 |

### レポート・モデル

| ファイル | 内容 |
|----------|------|
| `reports/mvp_report.md` | 自動生成 MVP レポート |
| `models/baseline_model.pkl` | 学習済みモデル |
| `models/model_metrics.json` | 評価指標 |
| `models/feature_importance.csv` | 特徴量重要度 |

---

## 8. 構造グループの分類ルール

| グループ | 条件 |
|----------|------|
| `G4_like` | G4Hunter風スコア ≥ 1.0 **または** G4 motif (`G3+N1-7G3+N1-7G3+N1-7G3+`) が存在 |
| `stem_loop_like` | paired base fraction > 0.3 **かつ** hairpin loop ≥ 1 |
| `G4_and_stem_loop_like` | 上記両方を満たす |
| `other_or_uncertain` | いずれも満たさない |

### G4 判定の詳細
- **G4Hunter 風スコア**: 連続 G に正スコア、連続 C に負スコアを付与し、  
  スライディングウィンドウ（デフォルト 25 bp）で平均し、最大絶対値をスコアとする
- **G4 motif 検出**: 正規表現 `G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}`

### ステムループ判定の詳細
- ViennaRNA (`RNA.fold`) が利用可能なら使用し、MFE 構造を取得
- 利用不可の場合、簡易 Nussinov アルゴリズムで相補対を推定
- dot-bracket から paired fraction / stem 数 / hairpin 数を抽出

---

## 9. タンパク質グループの分類ルール

### 電荷グループ (`charge_group`)

| グループ | 条件 |
|----------|------|
| `positive_charge_like` | estimated_charge_pH7.4 ≥ 5 |
| `negative_charge_like` | estimated_charge_pH7.4 ≤ -5 |
| `neutral_charge_like` | その間 |

### 疎水性グループ (`hydrophobicity_group`)

| グループ | 条件 |
|----------|------|
| `hydrophobic_like` | GRAVY ≥ 0.2 |
| `hydrophilic_like` | GRAVY ≤ -0.5 |
| `intermediate` | その間 |

### pI グループ (`pI_group`)

| グループ | 条件 |
|----------|------|
| `acidic_pI` | pI < 6 |
| `neutral_pI` | 6 ≤ pI ≤ 8 |
| `basic_pI` | pI > 8 |

### 機能グループ (`functional_group`)

UniProt keyword / GO term / subcellular location / タンパク質名に基づく:

| グループ | キーワード例 |
|----------|-------------|
| `membrane_protein` | membrane, transmembrane |
| `secreted_or_extracellular` | secreted, extracellular |
| `receptor` | receptor |
| `growth_factor_or_cytokine` | cytokine, growth factor, interleukin |
| `enzyme` | enzyme, kinase, protease |
| `nucleic_acid_binding_protein` | dna-binding, rna-binding |
| `cytosolic_or_intracellular` | cytoplasm, nucleus |
| `other_or_unknown` | 上記に該当しない |

---

## 10. ベースラインモデルについて

> **この結果は探索的解析であり、実験的検証が必要です。**

- **非結合データを使っていない**ため、「結合する / しない」の予測ではない
- 正例データのみからタンパク質グループの傾向を学習する探索的モデル
- Logistic Regression / Random Forest / Gradient Boosting を比較
- 評価: accuracy, balanced accuracy, macro F1, cross-validation, confusion matrix
- Feature importance を CSV 出力

---

## 11. プロジェクト構成

```
aptamer-protein-model/
  README.md
  requirements.txt
  run_pipeline.py
  .gitignore
  data/
    raw/
      aptamers_sample.csv
    processed/
  notebooks/
  src/
    __init__.py
    data/
      __init__.py
      load_aptamer_data.py
      normalize_targets.py
      fetch_uniprot.py
    features/
      __init__.py
      aptamer_features.py
      structure_grouping.py
      protein_features.py
      protein_grouping.py
    analysis/
      __init__.py
      association_analysis.py
      visualization.py
    models/
      __init__.py
      train_baseline.py
      evaluate.py
    utils/
      __init__.py
      io.py
      logging.py
  reports/
    figures/
    mvp_report.md
  models/
  tests/
    __init__.py
    test_aptamer_features.py
    test_protein_features.py
```

---

## 12. 今後の拡張案

1. AptaDB / UTexas Aptamer Database / Apta-Index などの大規模データ統合
2. UniProt マッピング精度改善（fuzzy matching等）
3. ViennaRNA / NUPACK ベースの高精度二次構造推定
4. Negative / decoy ペア設計とバイナリ分類
5. Multi-label / hierarchical classification
6. SHAP による説明可能性追加
7. 未知アプタマー配列入力に対する推定 CLI / GUI
8. Web アプリケーション化

---

## ライセンス

（プロジェクトに応じて設定してください）

---

## 引用

アプタマーデータベース:
- AptaDB: https://rnajournal.cshlp.org/content/30/3/189.full
- UTexas Aptamer Database: https://sites.utexas.edu/aptamerdatabase/
- Apta-Index

G4 予測:
- G4Hunter: https://pmc.ncbi.nlm.nih.gov/articles/PMC6748775/

二次構造予測:
- ViennaRNA: https://github.com/ViennaRNA/ViennaRNA

タンパク質物性:
- UniProt: https://www.uniprot.org/
- ProtParam: https://web.expasy.org/protparam/
