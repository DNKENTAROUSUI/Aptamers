# Aptamer-Protein MVP Report

## 1. 目的
既知のタンパク質結合アプタマー正例から、アプタマー配列・構造様特徴とタンパク質物性・機能グループの関係を探索的に解析する。

## 2. 使用データ
- アプタマー数: 55
- タンパク質標的数: 36
- DNA/RNA比率: {'DNA': 35, 'RNA': 20}
- Kd値あり: 51
- データベース別件数: {'AptaDB': 32, 'UTexas': 11, 'Ribocentre': 7, 'Apta-Index': 5}

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
| charge_group | 4.185 | 0.6517 | 6 |
| hydrophobicity_group | 7.643 | 0.2654 | 6 |
| pI_group | 5.153 | 0.5243 | 6 |
| functional_group | 28.44 | 0.004773 | 12 |

Fisher exact test summary:

| protein_group_column | aptamer_structure_group | protein_group | count_in_cell | odds_ratio | p_value | p_value_fdr_bh |
| --- | --- | --- | --- | --- | --- | --- |
| functional_group | stem_loop_like | nucleic_acid_binding_protein | 0 | 0 | 0.01886 | 0.198 |
| functional_group | other_or_uncertain | nucleic_acid_binding_protein | 2 | inf | 0.0101 | 0.198 |
| hydrophobicity_group | stem_loop_like | intermediate | 19 | inf | 0.0316 | 0.2844 |
| hydrophobicity_group | other_or_uncertain | intermediate | 0 | 0 | 0.06336 | 0.2851 |
| functional_group | other_or_uncertain | receptor | 2 | 11.75 | 0.05464 | 0.3825 |
| functional_group | stem_loop_like | membrane_protein | 15 | inf | 0.09095 | 0.406 |
| functional_group | stem_loop_like | receptor | 2 | 0.1333 | 0.09667 | 0.406 |
| hydrophobicity_group | other_or_uncertain | hydrophilic_like | 3 | 3.955 | 0.1656 | 0.4967 |
| functional_group | other_or_uncertain | membrane_protein | 0 | 0 | 0.1734 | 0.5203 |
| functional_group | G4_and_stem_loop_like | other_or_unknown | 2 | inf | 0.1704 | 0.5203 |
| hydrophobicity_group | stem_loop_like | hydrophilic_like | 11 | 0.3929 | 0.3564 | 0.5346 |
| hydrophobicity_group | stem_loop_like | hydrophobic_like | 9 | 0.3 | 0.3188 | 0.5346 |
| hydrophobicity_group | G4_and_stem_loop_like | hydrophobic_like | 1 | inf | 0.2667 | 0.5346 |
| hydrophobicity_group | other_or_uncertain | hydrophobic_like | 2 | 2 | 0.5981 | 0.769 |
| pI_group | other_or_uncertain | neutral_pI | 0 | 0 | 0.3131 | 0.8 |
| pI_group | stem_loop_like | basic_pI | 11 | 0.3929 | 0.3564 | 0.8 |
| pI_group | G4_and_stem_loop_like | acidic_pI | 1 | inf | 0.4444 | 0.8 |
| pI_group | other_or_uncertain | basic_pI | 3 | 3.955 | 0.1656 | 0.8 |
| pI_group | stem_loop_like | neutral_pI | 11 | inf | 0.3111 | 0.8 |
| charge_group | stem_loop_like | positive_charge_like | 11 | 0.3929 | 0.3564 | 1 |

自然言語サマリー:
- FDR補正後p<0.05で明確に過剰出現する構造群とタンパク質群の組み合わせは見つかりませんでした。

## 10. ベースラインモデル結果
目的変数 `hydrophobicity_group` に対して `gradient_boosting` を最良モデルとして保存した。

| model | cv_accuracy_mean | cv_balanced_accuracy_mean | cv_macro_f1_mean | accuracy | balanced_accuracy | macro_f1 |
| --- | --- | --- | --- | --- | --- | --- |
| logistic_regression | 0.5476190476190477 | 0.5777777777777777 | 0.48698412698412696 | 0.5833333333333334 | 0.5888888888888889 | 0.6151515151515151 |
| random_forest | 0.780952380952381 | 0.8 | 0.7787301587301588 | 0.75 | 0.7666666666666666 | 0.7428571428571429 |
| gradient_boosting | 0.9428571428571428 | 0.9333333333333332 | 0.9166666666666667 | 1.0 | 1.0 | 1.0 |

この結果は探索的解析であり、実験的検証が必要である。

## 11. 限界
この結果は既知の正例データに基づく探索的解析であり、非結合を直接予測するものではない。外部DB由来データには選択バイアス、表記揺れ、Kd条件差、配列切り出し差がある。二次構造は予測であり、実験的検証が必要である。

## 12. 次にやるべきこと
データソースを拡充し、標的名の正規化とUniProt accession対応を改善する。測定条件を標準化し、負例または候補集合を設計したうえで、独立テストセットによる予測性能評価を行う。
