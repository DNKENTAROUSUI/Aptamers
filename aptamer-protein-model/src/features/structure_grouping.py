from __future__ import annotations

G4_MIN_G_FRACTION = 0.30
G4_MIN_GRUN_COUNT = 3
G4_MIN_HUNTER_SCORE = 1.2
G4_MOTIF_REGEX = r"G{3,}[ACGTU]{1,7}G{3,}[ACGTU]{1,7}G{3,}[ACGTU]{1,7}G{3,}"

STEM_LOOP_MIN_PAIRED_FRACTION = 0.20
STEM_LOOP_MIN_STEMS = 1
STEM_LOOP_MIN_HAIRPINS = 1


def classify_g4_like(g_fraction: float, g_run_count: int, g4hunter_score: float, has_g4_motif: bool) -> str:
    is_g4 = bool(has_g4_motif) or (
        g_fraction >= G4_MIN_G_FRACTION
        and g_run_count >= G4_MIN_GRUN_COUNT
        and g4hunter_score >= G4_MIN_HUNTER_SCORE
    )
    return "G4_like" if is_g4 else "non_G4_like"


def classify_stem_loop_like(paired_fraction: float, stem_count: int, hairpin_loop_count: int) -> str:
    is_stem_loop = (
        paired_fraction >= STEM_LOOP_MIN_PAIRED_FRACTION
        and stem_count >= STEM_LOOP_MIN_STEMS
        and hairpin_loop_count >= STEM_LOOP_MIN_HAIRPINS
    )
    return "stem_loop_like" if is_stem_loop else "non_stem_loop_like"


def classify_structure_group(g4_group: str, stem_loop_group: str) -> str:
    g4 = g4_group == "G4_like"
    stem = stem_loop_group == "stem_loop_like"
    if g4 and stem:
        return "G4_and_stem_loop_like"
    if g4:
        return "G4_like"
    if stem:
        return "stem_loop_like"
    return "other_or_uncertain"
