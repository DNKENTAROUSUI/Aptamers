"""
アプタマー特徴量のユニットテスト
"""

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd
from src.features.aptamer_features import (
    compute_basic_features,
    compute_g4_features,
    compute_kmer_features,
    compute_stem_loop_features,
    classify_structure_group,
    compute_all_features,
)


def test_basic_features_thrombin_aptamer():
    seq = "GGTTGGTGTGGTTGG"
    feat = compute_basic_features(seq, "DNA")
    assert feat["sequence_length"] == 15
    assert feat["gc_content"] >= 0
    assert feat["max_consecutive_G"] >= 2
    assert feat["nucleic_acid_type_detected"] == "DNA"


def test_basic_features_rna():
    seq = "CGGAAUCAGUGAAUGCUUAUACAUCCG"
    feat = compute_basic_features(seq, "RNA")
    assert feat["nucleic_acid_type_detected"] == "RNA"
    assert feat["frac_U"] > 0
    assert feat["frac_T"] == 0


def test_basic_features_empty():
    feat = compute_basic_features("", "")
    assert feat["sequence_length"] == 0


def test_kmer_features():
    seq = "GGTTGGTGTGGTTGG"
    feat = compute_kmer_features(seq, k_values=(1, 2))
    assert "kmer_G" in feat
    assert "kmer_GG" in feat
    assert feat["kmer_G"] > 0


def test_g4_features_g_rich():
    seq = "GGGTTAGGGTTAGGGTTAGGG"
    feat = compute_g4_features(seq)
    assert feat["g_richness"] > 0.4
    assert feat["n_g_runs"] >= 4
    assert feat["has_g4_motif"] == 1
    assert feat["is_g4_like"] == 1


def test_g4_features_non_g_rich():
    seq = "AATTCCAATTCCAATTCC"
    feat = compute_g4_features(seq)
    assert feat["g_richness"] == 0.0
    assert feat["is_g4_like"] == 0


def test_stem_loop_features():
    # Self-complementary region should form pairs
    seq = "GGGGAAAACCCC"
    feat = compute_stem_loop_features(seq, "DNA")
    assert "paired_base_fraction" in feat
    assert "mfe" in feat


def test_classify_structure_group():
    assert classify_structure_group(True, False) == "G4_like"
    assert classify_structure_group(False, True) == "stem_loop_like"
    assert classify_structure_group(True, True) == "G4_and_stem_loop_like"
    assert classify_structure_group(False, False) == "other_or_uncertain"


def test_compute_all_features():
    df = pd.DataFrame({
        "aptamer_id": ["A1", "A2"],
        "sequence": ["GGGTTAGGGTTAGGGTTAGGG", "AATTCCGGAATTCCGG"],
        "nucleic_acid_type": ["DNA", "DNA"],
        "target_name": ["Thrombin", "Lysozyme"],
    })
    out = compute_all_features(df)
    assert len(out) == 2
    assert "structure_group" in out.columns
    assert "g4_hunter_score" in out.columns
    assert "mfe" in out.columns
    assert "gc_content" in out.columns


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
