# Aptamer-Protein MVP Report

## 1. 目的
既知のタンパク質結合アプタマー正例から、アプタマー配列・構造様特徴とタンパク質物性・機能グループの関係を探索的に解析する。

## 2. 使用データ
- アプタマー数: 162
- タンパク質標的数: 36
- DNA/RNA比率: {'DNA': 102, 'RNA': 60}
- Kd値あり: 50
- データベース別件数: {'exploratory_variant': 108, 'AptaDB': 32, 'UTexas': 11, 'Ribocentre': 7, 'Apta-Index': 4}
- 候補ステータス別件数: {'in_silico_variant': 108, 'known_or_reported': 54}

## 3. データ取得方法
MVPではローカルCSVを優先する。タンパク質配列はローカルFASTA、または `--fetch-uniprot` 指定時にUniProt REST APIから取得する。`in_silico_variant` は既知結合の正例ではなく、探索用の派生候補として扱う。

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
| charge_group | 5.304 | 0.5054 | 6 |
| hydrophobicity_group | 23.73 | 0.0005847 | 6 |
| pI_group | 8.897 | 0.1795 | 6 |
| functional_group | 85.1 | 4.36e-13 | 12 |

Fisher exact test summary:

| protein_group_column | aptamer_structure_group | protein_group | count_in_cell | odds_ratio | p_value | p_value_fdr_bh |
| --- | --- | --- | --- | --- | --- | --- |
| functional_group | other_or_uncertain | nucleic_acid_binding_protein | 6 | inf | 7.506e-08 | 1.576e-06 |
| functional_group | stem_loop_like | nucleic_acid_binding_protein | 0 | 0 | 8.121e-07 | 8.527e-06 |
| hydrophobicity_group | stem_loop_like | intermediate | 57 | inf | 0.0002746 | 0.002472 |
| hydrophobicity_group | other_or_uncertain | hydrophilic_like | 9 | 12 | 0.0005555 | 0.0025 |
| hydrophobicity_group | other_or_uncertain | intermediate | 0 | 0 | 0.002404 | 0.007213 |
| functional_group | stem_loop_like | membrane_protein | 45 | inf | 0.003674 | 0.02572 |
| hydrophobicity_group | G4_and_stem_loop_like | hydrophobic_like | 3 | inf | 0.01456 | 0.02621 |
| hydrophobicity_group | stem_loop_like | hydrophilic_like | 33 | 0.2157 | 0.01203 | 0.02621 |
| functional_group | G4_and_stem_loop_like | other_or_unknown | 5 | inf | 0.01286 | 0.06754 |
| functional_group | other_or_uncertain | membrane_protein | 0 | 0 | 0.02035 | 0.08546 |
| functional_group | other_or_uncertain | other_or_unknown | 2 | 0.2225 | 0.04381 | 0.1533 |
| pI_group | stem_loop_like | neutral_pI | 33 | inf | 0.02071 | 0.1864 |
| pI_group | other_or_uncertain | neutral_pI | 0 | 0 | 0.06464 | 0.274 |
| pI_group | G4_and_stem_loop_like | acidic_pI | 3 | inf | 0.09134 | 0.274 |
| pI_group | stem_loop_like | acidic_pI | 51 | 0.4229 | 0.1621 | 0.3648 |
| hydrophobicity_group | G4_and_stem_loop_like | intermediate | 0 | 0 | 0.2583 | 0.3875 |
| hydrophobicity_group | stem_loop_like | hydrophobic_like | 28 | 0.56 | 0.3378 | 0.4343 |
| functional_group | other_or_uncertain | receptor | 2 | 3.688 | 0.1553 | 0.4659 |
| pI_group | other_or_uncertain | basic_pI | 5 | 2.132 | 0.3 | 0.54 |
| pI_group | G4_and_stem_loop_like | basic_pI | 0 | 0 | 0.5547 | 0.5727 |

自然言語サマリー:
- other_or_uncertain は nucleic_acid_binding_protein (functional_group) に多い傾向があります (odds ratio=inf, FDR p=1.58e-06)。
- stem_loop_like は intermediate (hydrophobicity_group) に多い傾向があります (odds ratio=inf, FDR p=0.00247)。
- other_or_uncertain は hydrophilic_like (hydrophobicity_group) に多い傾向があります (odds ratio=12.00, FDR p=0.0025)。
- stem_loop_like は membrane_protein (functional_group) に多い傾向があります (odds ratio=inf, FDR p=0.0257)。
- G4_and_stem_loop_like は hydrophobic_like (hydrophobicity_group) に多い傾向があります (odds ratio=inf, FDR p=0.0262)。

## 10. ベースラインモデル結果
目的変数 `hydrophobicity_group` に対して `gradient_boosting` を最良モデルとして保存した。

| model | cv_accuracy_mean | cv_balanced_accuracy_mean | cv_macro_f1_mean | accuracy | balanced_accuracy | macro_f1 |
| --- | --- | --- | --- | --- | --- | --- |
| logistic_regression | 0.5476190476190477 | 0.5777777777777777 | 0.48698412698412696 | 0.5833333333333334 | 0.5888888888888889 | 0.6151515151515151 |
| random_forest | 0.780952380952381 | 0.8 | 0.7787301587301588 | 0.75 | 0.7666666666666666 | 0.7428571428571429 |
| gradient_boosting | 0.9428571428571428 | 0.9333333333333332 | 0.9166666666666667 | 1.0 | 1.0 | 1.0 |

この結果は探索的解析であり、実験的検証が必要である。

## 11. 限界
この結果は既知の正例データと探索用派生候補に基づく解析であり、非結合を直接予測するものではない。`in_silico_variant` は実験的に結合確認された配列ではない。外部DB由来データには選択バイアス、表記揺れ、Kd条件差、配列切り出し差がある。二次構造は予測であり、実験的検証が必要である。

## 12. 次にやるべきこと
データソースを拡充し、標的名の正規化とUniProt accession対応を改善する。測定条件を標準化し、負例または候補集合を設計したうえで、独立テストセットによる予測性能評価を行う。
