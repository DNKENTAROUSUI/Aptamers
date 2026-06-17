from __future__ import annotations

CHARGE_POSITIVE_THRESHOLD = 1.0
CHARGE_NEGATIVE_THRESHOLD = -1.0
HYDROPHOBIC_GRAVY_THRESHOLD = 0.0
HYDROPHILIC_GRAVY_THRESHOLD = -0.4
ACIDIC_PI_THRESHOLD = 6.0
BASIC_PI_THRESHOLD = 8.0

FUNCTION_KEYWORDS = {
    "membrane_protein": ("membrane", "transmembrane", "cell membrane"),
    "secreted_or_extracellular": ("secreted", "extracellular", "signal peptide"),
    "cytosolic_or_intracellular": ("cytoplasm", "cytosol", "intracellular"),
    "enzyme": ("enzyme", "kinase", "protease", "hydrolase", "transferase", "oxidoreductase"),
    "receptor": ("receptor",),
    "growth_factor_or_cytokine": ("growth factor", "cytokine", "interleukin", "chemokine"),
    "nucleic_acid_binding_protein": ("dna-binding", "rna-binding", "nucleic acid-binding", "transcription factor"),
}


def charge_group(charge_at_7_4: float | None) -> str:
    if charge_at_7_4 is None or charge_at_7_4 != charge_at_7_4:
        return "unknown"
    if charge_at_7_4 > CHARGE_POSITIVE_THRESHOLD:
        return "positive_charge_like"
    if charge_at_7_4 < CHARGE_NEGATIVE_THRESHOLD:
        return "negative_charge_like"
    return "neutral_charge_like"


def hydrophobicity_group(gravy: float | None) -> str:
    if gravy is None or gravy != gravy:
        return "unknown"
    if gravy > HYDROPHOBIC_GRAVY_THRESHOLD:
        return "hydrophobic_like"
    if gravy < HYDROPHILIC_GRAVY_THRESHOLD:
        return "hydrophilic_like"
    return "intermediate"


def pi_group(pi: float | None) -> str:
    if pi is None or pi != pi:
        return "unknown"
    if pi < ACIDIC_PI_THRESHOLD:
        return "acidic_pI"
    if pi > BASIC_PI_THRESHOLD:
        return "basic_pI"
    return "neutral_pI"


def functional_group(annotation_text: str) -> str:
    text = annotation_text.lower()
    for group, keywords in FUNCTION_KEYWORDS.items():
        if any(keyword in text for keyword in keywords):
            return group
    return "other_or_unknown"
