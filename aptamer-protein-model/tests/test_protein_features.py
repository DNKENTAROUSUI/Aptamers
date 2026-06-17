from src.features.protein_features import compute_physical_features
from src.features.protein_grouping import charge_group, hydrophobicity_group, pi_group


def test_compute_physical_features_basic():
    features = compute_physical_features("MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE")
    assert features["protein_sequence_length"] == 40
    assert features["molecular_weight"] > 0
    assert "GRAVY" in features


def test_grouping_rules():
    assert charge_group(3.0) == "positive_charge_like"
    assert charge_group(-3.0) == "negative_charge_like"
    assert hydrophobicity_group(0.5) == "hydrophobic_like"
    assert pi_group(9.0) == "basic_pI"
