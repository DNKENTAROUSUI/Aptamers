# Aptamer-Protein MVP Report

## 1. 目的
既知のタンパク質結合アプタマー正例から、アプタマー配列・構造様特徴とタンパク質物性・機能グループの関係を探索的に解析する。

## 2. 使用データ
- アプタマー数: 5
- タンパク質標的数: 5
- DNA/RNA比率: {'DNA': 4, 'RNA': 1}
- Kd値あり: 5
- データベース別件数: {'example': 5}

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

| protein_group_column | chi2 | p_value | dof |
| --- | --- | --- | --- |
| charge_group | 5 | 0.08208 | 2 |
| hydrophobicity_group | nan | nan | nan |
| pI_group | 5 | 0.08208 | 2 |
| functional_group | nan | nan | nan |

Fisher exact test summary:

| protein_group_column | aptamer_structure_group | protein_group | count_in_cell | odds_ratio | p_value | p_value_fdr_bh |
| --- | --- | --- | --- | --- | --- | --- |
| charge_group | other_or_uncertain | negative_charge_like | 2 | inf | 0.1 | 0.3 |
| charge_group | stem_loop_like | negative_charge_like | 0 | 0 | 0.1 | 0.3 |
| pI_group | stem_loop_like | acidic_pI | 0 | 0 | 0.1 | 0.3 |
| pI_group | other_or_uncertain | acidic_pI | 2 | inf | 0.1 | 0.3 |
| pI_group | stem_loop_like | neutral_pI | 2 | inf | 0.4 | 0.6 |
| pI_group | other_or_uncertain | neutral_pI | 0 | 0 | 0.4 | 0.6 |
| charge_group | other_or_uncertain | positive_charge_like | 0 | 0 | 0.4 | 0.6 |
| charge_group | stem_loop_like | positive_charge_like | 2 | inf | 0.4 | 0.6 |
| hydrophobicity_group | stem_loop_like | hydrophilic_like | 3 | nan | 1 | 1 |
| hydrophobicity_group | other_or_uncertain | hydrophilic_like | 2 | nan | 1 | 1 |
| charge_group | stem_loop_like | neutral_charge_like | 1 | inf | 1 | 1 |
| charge_group | other_or_uncertain | neutral_charge_like | 0 | 0 | 1 | 1 |
| pI_group | other_or_uncertain | basic_pI | 0 | 0 | 1 | 1 |
| pI_group | stem_loop_like | basic_pI | 1 | inf | 1 | 1 |
| functional_group | other_or_uncertain | other_or_unknown | 2 | nan | 1 | 1 |
| functional_group | stem_loop_like | other_or_unknown | 3 | nan | 1 | 1 |

自然言語サマリー:
- FDR補正後p<0.05で明確に過剰出現する構造群とタンパク質群の組み合わせは見つかりませんでした。

## 10. ベースラインモデル結果
モデル学習は未実行、またはデータ数不足のためスキップされました。

## 11. 限界
この結果は既知の正例データに基づく探索的解析であり、非結合を直接予測するものではない。外部DB由来データには選択バイアス、表記揺れ、Kd条件差、配列切り出し差がある。二次構造は予測であり、実験的検証が必要である。

## 12. 次にやるべきこと
データソースを拡充し、標的名の正規化とUniProt accession対応を改善する。測定条件を標準化し、負例または候補集合を設計したうえで、独立テストセットによる予測性能評価を行う。
