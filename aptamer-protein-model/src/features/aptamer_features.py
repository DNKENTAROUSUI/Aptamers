from __future__ import annotations

import argparse
import itertools
import math
import re
import subprocess
from collections import Counter

import numpy as np
import pandas as pd

from src.features.structure_grouping import (
    G4_MOTIF_REGEX,
    classify_g4_like,
    classify_stem_loop_like,
    classify_structure_group,
)
from src.utils.io import ensure_parent, read_csv, write_csv
from src.utils.logging import get_logger

LOGGER = get_logger(__name__)
BASES = "ACGTU"
PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T", "U"}


def clean_sequence(seq: object) -> str:
    if pd.isna(seq):
        return ""
    return re.sub(r"[^ACGTU]", "", str(seq).upper().replace(" ", "").replace("\n", ""))


def max_run(seq: str, base: str) -> int:
    runs = re.findall(f"{base}+", seq)
    return max((len(run) for run in runs), default=0)


def kmer_frequencies(seq: str, k: int) -> dict[str, float]:
    alphabet = "ACGU" if "U" in seq and "T" not in seq else "ACGT"
    keys = ["".join(kmer) for kmer in itertools.product(alphabet, repeat=k)]
    total = max(len(seq) - k + 1, 0)
    counts = Counter(seq[i : i + k] for i in range(total))
    return {f"{k}mer_{key}": counts.get(key, 0) / total if total else 0.0 for key in keys}


def g4hunter_like_score(seq: str, window: int = 25) -> float:
    if not seq:
        return np.nan
    scores = {"G": 1.0, "C": -1.0, "A": 0.0, "T": 0.0, "U": 0.0}
    values = [scores.get(base, 0.0) for base in seq]
    if len(values) <= window:
        return float(abs(np.mean(values)))
    return float(max(abs(np.mean(values[i : i + window])) for i in range(len(values) - window + 1)))


def g_run_loop_lengths(seq: str) -> tuple[int, float, int]:
    matches = list(re.finditer(r"G{2,}", seq))
    if len(matches) < 2:
        return len(matches), np.nan, 0
    loops = [matches[i + 1].start() - matches[i].end() for i in range(len(matches) - 1)]
    return len(matches), float(np.mean(loops)), int(max(loops))


def can_pair(a: str, b: str) -> bool:
    return (a, b) in {("A", "T"), ("T", "A"), ("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}


def nussinov_dot_bracket(seq: str, min_loop: int = 3) -> str:
    n = len(seq)
    if n == 0:
        return ""
    dp = [[0] * n for _ in range(n)]
    trace: list[tuple[int, int]] = []
    for length in range(1, n):
        for i in range(n - length):
            j = i + length
            best = dp[i + 1][j]
            best = max(best, dp[i][j - 1])
            if j - i > min_loop and can_pair(seq[i], seq[j]):
                best = max(best, dp[i + 1][j - 1] + 1)
            for k in range(i + 1, j):
                best = max(best, dp[i][k] + dp[k + 1][j])
            dp[i][j] = best

    pairs: list[tuple[int, int]] = []

    def backtrack(i: int, j: int) -> None:
        if i >= j:
            return
        if dp[i][j] == dp[i + 1][j]:
            backtrack(i + 1, j)
        elif dp[i][j] == dp[i][j - 1]:
            backtrack(i, j - 1)
        elif j - i > min_loop and can_pair(seq[i], seq[j]) and dp[i][j] == dp[i + 1][j - 1] + 1:
            pairs.append((i, j))
            backtrack(i + 1, j - 1)
        else:
            for k in range(i + 1, j):
                if dp[i][j] == dp[i][k] + dp[k + 1][j]:
                    backtrack(i, k)
                    backtrack(k + 1, j)
                    break

    backtrack(0, n - 1)
    chars = ["."] * n
    for i, j in pairs:
        chars[i] = "("
        chars[j] = ")"
    return "".join(chars)


def run_rnafold(seq: str) -> tuple[str, float] | None:
    try:
        completed = subprocess.run(
            ["RNAfold", "--noPS"],
            input=f"{seq}\n",
            text=True,
            capture_output=True,
            check=True,
            timeout=10,
        )
    except (FileNotFoundError, subprocess.SubprocessError):
        return None
    lines = completed.stdout.strip().splitlines()
    if len(lines) < 2:
        return None
    match = re.search(r"([().]+)\s+\(([-0-9.]+)\)", lines[1])
    if not match:
        return None
    return match.group(1), float(match.group(2))


def parse_dot_bracket(dot: str) -> dict[str, float | int]:
    paired = dot.count("(") + dot.count(")")
    paired_fraction = paired / len(dot) if dot else 0.0
    stems = re.findall(r"\(+", dot)
    loops = []
    for match in re.finditer(r"\.+", dot):
        left = dot[: match.start()]
        right = dot[match.end() :]
        if left.endswith("(") and right.startswith(")"):
            loops.append(match.group(0))
    all_unpaired = re.findall(r"\.+", dot)
    return {
        "paired_base_fraction": paired_fraction,
        "stem_count": len(stems),
        "hairpin_loop_count": len(loops),
        "max_stem_length": max((len(s) for s in stems), default=0),
        "max_loop_length": max((len(loop) for loop in all_unpaired), default=0),
    }


def infer_nucleic_acid_type(seq: str, provided: object = None) -> str:
    if isinstance(provided, str) and provided.strip():
        value = provided.strip().upper()
        if value in {"DNA", "RNA"}:
            return value
    if "U" in seq and "T" not in seq:
        return "RNA"
    if "T" in seq and "U" not in seq:
        return "DNA"
    return "unknown"


def compute_aptamer_features(row: pd.Series) -> dict[str, object]:
    seq = clean_sequence(row.get("sequence", ""))
    length = len(seq)
    counts = Counter(seq)
    freqs = {f"{base}_fraction": counts.get(base, 0) / length if length else 0.0 for base in BASES}
    gc = (counts.get("G", 0) + counts.get("C", 0)) / length if length else 0.0
    purine = sum(counts.get(b, 0) for b in PURINES) / length if length else 0.0
    pyrimidine = sum(counts.get(b, 0) for b in PYRIMIDINES) / length if length else 0.0
    grun_count, mean_loop, max_loop = g_run_loop_lengths(seq)
    g4_score = g4hunter_like_score(seq)
    g_fraction = freqs["G_fraction"]
    c_fraction = freqs["C_fraction"]
    g_skewness = (g_fraction - c_fraction) / (g_fraction + c_fraction) if (g_fraction + c_fraction) else 0.0
    has_g4_motif = bool(re.search(G4_MOTIF_REGEX, seq))

    folded = run_rnafold(seq.replace("T", "U"))
    if folded:
        dot, mfe = folded
    else:
        dot = nussinov_dot_bracket(seq)
        mfe = -1.0 * dot.count("(")
    structure = parse_dot_bracket(dot)
    g4_group = classify_g4_like(g_fraction, grun_count, g4_score, has_g4_motif)
    stem_group = classify_stem_loop_like(
        structure["paired_base_fraction"], structure["stem_count"], structure["hairpin_loop_count"]
    )
    result: dict[str, object] = {
        "aptamer_id": row.get("aptamer_id"),
        "aptamer_name": row.get("aptamer_name"),
        "sequence": seq,
        "nucleic_acid_type": infer_nucleic_acid_type(seq, row.get("nucleic_acid_type")),
        "sequence_length": length,
        **freqs,
        "GC_content": gc,
        "G_content": g_fraction,
        "C_content": c_fraction,
        "purine_fraction": purine,
        "pyrimidine_fraction": pyrimidine,
        "max_consecutive_G": max_run(seq, "G"),
        "max_consecutive_C": max_run(seq, "C"),
        "G_richness": g_fraction,
        "G_skewness": g_skewness,
        "G4Hunter_like_score": g4_score,
        "has_G_quadruplex_motif": has_g4_motif,
        "G_run_count": grun_count,
        "mean_G_run_loop_length": mean_loop,
        "max_G_run_loop_length": max_loop,
        "MFE": mfe,
        "predicted_dot_bracket": dot,
        **structure,
        "g4_group": g4_group,
        "stem_loop_group": stem_group,
        "aptamer_structure_group": classify_structure_group(g4_group, stem_group),
    }
    for k in (1, 2, 3):
        result.update(kmer_frequencies(seq, k))
    return result


def build_aptamer_feature_table(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([compute_aptamer_features(row) for _, row in df.iterrows()])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/processed/aptamer_protein_pairs_clean.csv")
    parser.add_argument("--output", default="data/processed/aptamer_features.csv")
    args = parser.parse_args()
    df = read_csv(args.input)
    features = build_aptamer_feature_table(df)
    write_csv(features, args.output)
    LOGGER.info("Wrote %s rows to %s", len(features), args.output)


if __name__ == "__main__":
    main()
