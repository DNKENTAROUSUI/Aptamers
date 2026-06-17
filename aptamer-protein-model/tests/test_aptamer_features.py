import pandas as pd

from src.features.aptamer_features import compute_aptamer_features


def test_compute_aptamer_features_g4_like():
    row = pd.Series(
        {
            "aptamer_id": "x",
            "aptamer_name": "g-rich",
            "sequence": "GGGAGGGTTGGGACGGG",
            "nucleic_acid_type": "DNA",
        }
    )
    features = compute_aptamer_features(row)
    assert features["sequence_length"] == 17
    assert features["has_G_quadruplex_motif"] is True
    assert features["g4_group"] == "G4_like"
    assert "aptamer_structure_group" in features


def test_kmer_features_present():
    features = compute_aptamer_features(pd.Series({"sequence": "ACGTACGT", "nucleic_acid_type": "DNA"}))
    assert "1mer_A" in features
    assert "2mer_AC" in features
    assert "3mer_ACG" in features
