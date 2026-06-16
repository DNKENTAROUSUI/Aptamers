# MVP Report: Aptamer-Protein Structure Association Analysis

> このレポートは `python -m src.analysis.association_analysis` 実行時に自動生成で上書きされます。

## 1. 目的
既知のタンパク質結合アプタマーの配列特徴と簡易構造特徴を抽出し、
タンパク質の物性・機能グループとの関係を探索的に解析する。

## 2. 使用データ
（自動生成で上書き）

## 3. データ取得方法
- ローカルCSVを入力
- target_type == protein を対象
- 必要に応じて UniProt API で補完

## 4. アプタマー特徴量
- 配列長、塩基組成、GC含量、k-mer
- G4関連特徴量
- ステムループ関連特徴量

## 5. タンパク質特徴量
- 長さ、分子量、pI、推定電荷、GRAVY
- アミノ酸組成群割合

## 6. 構造グループの分類ルール
- G4_like / stem_loop_like / G4_and_stem_loop_like / other_or_uncertain

## 7. タンパク質グループの分類ルール
- charge_group / hydrophobicity_group / pI_group / functional_group

## 8–12. （自動生成で上書き）

## 限界
- 非結合データなし
- 探索的解析であり、実験的検証が必要
