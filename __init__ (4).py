"""
アプタマー配列特徴量計算モジュール
======================================
"""

import argparse
import math
import os
import re
import sys
import warnings
from collections import Counter
from itertools import product
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from src.utils.logging import setup_logger
from src.utils.io import safe_read_csv, safe_write_csv

logger = setup_logger(__name__)

# ================================================================
# 定数・閾値
# ================================================================
G4_HUNTER_WINDOW = 25
G4_HUNTER_THRESHOLD = 1.0
G4_MOTIF_MIN_G_RUN = 3
G4_MOTIF_MAX_LOOP = 7
G4_MOTIF_MIN_LOOP = 1

STEM_LOOP_PAIRED_FRACTION_THRESHOLD = 0.3
STEM_LOOP_MIN_HAIRPIN = 1

G4_MOTIF_PATTERN = re.compile(
    r"[Gg]{3,}[^Gg]{1,7}[Gg]{3,}[^Gg]{1,7}[Gg]{3,}[^Gg]{1,7}[Gg]{3,}"
)

_COMPLEMENT_DNA = {"A": "T", "T": "A", "G": "C", "C": "G"}
_COMPLEMENT_RNA = {"A": "U", "U": "A", "G": "C", "C": "G"}


# ================================================================
# 基本特徴量
# ================================================================

def detect_nucleic_acid_type(sequence, annotated_type=""):
    if annotated_type and annotated_type.upper() in ("DNA", "RNA"):
        return annotated_type.upper()
    seq_upper = sequence.upper()
    if "U" in seq_upper and "T" not in seq_upper:
        return "RNA"
    elif "T" in seq_upper and "U" not in seq_upper:
        return "DNA"
    return "DNA"


def compute_basic_features(sequence, na_type=""):
    seq = sequence.upper().strip()
    length = len(seq)
    if length == 0:
        return _empty_basic_features()
    na = detect_nucleic_acid_type(seq, na_type)
    counts = Counter(seq)
    a_frac = counts.get("A", 0) / length
    c_frac = counts.get("C", 0) / length
    g_frac = counts.get("G", 0) / length
    t_frac = counts.get("T", 0) / length
    u_frac = counts.get("U", 0) / length
    gc_content = (counts.get("G", 0) + counts.get("C", 0)) / length
    purine_frac = (counts.get("A", 0) + counts.get("G", 0)) / length
    pyrimidine_frac = 1.0 - purine_frac
    max_g_run = _max_consecutive(seq, "G")
    max_c_run = _max_consecutive(seq, "C")
    return {
        "sequence_length": length,
        "nucleic_acid_type_detected": na,
        "frac_A": round(a_frac, 4), "frac_C": round(c_frac, 4),
        "frac_G": round(g_frac, 4), "frac_T": round(t_frac, 4),
        "frac_U": round(u_frac, 4),
        "gc_content": round(gc_content, 4),
        "g_content": round(g_frac, 4), "c_content": round(c_frac, 4),
        "purine_fraction": round(purine_frac, 4),
        "pyrimidine_fraction": round(pyrimidine_frac, 4),
        "max_consecutive_G": max_g_run, "max_consecutive_C": max_c_run,
    }


def _empty_basic_features():
    return {k: 0 for k in [
        "sequence_length", "frac_A", "frac_C", "frac_G", "frac_T", "frac_U",
        "gc_content", "g_content", "c_content", "purine_fraction",
        "pyrimidine_fraction", "max_consecutive_G", "max_consecutive_C",
    ]}


def _max_consecutive(seq, char):
    max_run = current = 0
    for c in seq:
        if c == char:
            current += 1
            max_run = max(max_run, current)
        else:
            current = 0
    return max_run


# ================================================================
# k-mer 特徴量
# ================================================================

def compute_kmer_features(sequence, k_values=(1, 2, 3)):
    seq = sequence.upper().strip()
    features = {}
    bases = "ACGTU"
    for k in k_values:
        if len(seq) < k:
            for kmer in product(bases, repeat=k):
                features[f"kmer_{''.join(kmer)}"] = 0.0
            continue
        total = len(seq) - k + 1
        kmer_counts = Counter(seq[i:i+k] for i in range(total))
        for kmer in product(bases, repeat=k):
            kmer_str = "".join(kmer)
            features[f"kmer_{kmer_str}"] = round(
                kmer_counts.get(kmer_str, 0) / total, 4
            ) if total > 0 else 0.0
    return features


# ================================================================
# G4 構造特徴量
# ================================================================

def compute_g4_features(sequence):
    seq = sequence.upper().strip()
    length = len(seq)
    if length == 0:
        return _empty_g4_features()
    g_count = seq.count("G")
    c_count = seq.count("C")
    g_richness = g_count / length
    g_skewness = (g_count - c_count) / (g_count + c_count) if (g_count + c_count) > 0 else 0.0
    g4_hunter_score = _g4_hunter_score(seq)
    has_g4_motif = bool(G4_MOTIF_PATTERN.search(seq))
    g_runs = re.findall(r"G{3,}", seq)
    n_g_runs = len(g_runs)
    avg_loop_length = _average_g4_loop_length(seq)
    is_g4_like = g4_hunter_score >= G4_HUNTER_THRESHOLD or has_g4_motif
    return {
        "g_richness": round(g_richness, 4),
        "g_skewness": round(g_skewness, 4),
        "g4_hunter_score": round(g4_hunter_score, 4),
        "has_g4_motif": int(has_g4_motif),
        "n_g_runs": n_g_runs,
        "avg_g4_loop_length": round(avg_loop_length, 2),
        "is_g4_like": int(is_g4_like),
    }


def _empty_g4_features():
    return {k: 0 for k in [
        "g_richness", "g_skewness", "g4_hunter_score", "has_g4_motif",
        "n_g_runs", "avg_g4_loop_length", "is_g4_like",
    ]}


def _g4_hunter_score(seq):
    if len(seq) == 0:
        return 0.0
    scores = np.zeros(len(seq))
    i = 0
    while i < len(seq):
        if seq[i] == "G":
            j = i
            while j < len(seq) and seq[j] == "G":
                j += 1
            run_len = j - i
            for k in range(i, j):
                scores[k] = min(run_len, 4)
            i = j
        elif seq[i] == "C":
            j = i
            while j < len(seq) and seq[j] == "C":
                j += 1
            run_len = j - i
            for k in range(i, j):
                scores[k] = -min(run_len, 4)
            i = j
        else:
            i += 1
    window = min(G4_HUNTER_WINDOW, len(seq))
    if window == 0:
        return 0.0
    max_abs_score = 0.0
    for start in range(len(seq) - window + 1):
        window_mean = np.mean(scores[start:start + window])
        max_abs_score = max(max_abs_score, abs(window_mean))
    return max_abs_score


def _average_g4_loop_length(seq):
    g_run_positions = [(m.start(), m.end()) for m in re.finditer(r"G{3,}", seq)]
    if len(g_run_positions) < 2:
        return 0.0
    loop_lengths = []
    for i in range(len(g_run_positions) - 1):
        loop_len = g_run_positions[i + 1][0] - g_run_positions[i][1]
        loop_lengths.append(loop_len)
    return np.mean(loop_lengths) if loop_lengths else 0.0


# ================================================================
# ステムループ構造特徴量
# ================================================================

def compute_stem_loop_features(sequence, na_type="DNA"):
    seq = sequence.upper().strip()
    length = len(seq)
    if length == 0:
        return _empty_stem_loop_features()
    vienna_result = _try_vienna_fold(seq)
    if vienna_result is not None:
        return vienna_result
    return _simplified_structure_features(seq, na_type)


def _try_vienna_fold(seq):
    try:
        import RNA
        rna_seq = seq.replace("T", "U")
        structure, mfe = RNA.fold(rna_seq)
        return _parse_dot_bracket(seq, structure, mfe)
    except ImportError:
        return None
    except Exception as e:
        warnings.warn(f"ViennaRNA エラー: {e}")
        return None


def _parse_dot_bracket(seq, structure, mfe):
    length = len(seq)
    n_paired = structure.count("(") + structure.count(")")
    paired_fraction = n_paired / length if length > 0 else 0.0
    stems, hairpins = _count_stems_and_hairpins(structure)
    max_stem, max_loop = _max_stem_and_loop(structure)
    is_stem_loop = (
        paired_fraction >= STEM_LOOP_PAIRED_FRACTION_THRESHOLD
        and hairpins >= STEM_LOOP_MIN_HAIRPIN
    )
    return {
        "mfe": round(mfe, 2),
        "paired_base_fraction": round(paired_fraction, 4),
        "n_stems": stems, "n_hairpin_loops": hairpins,
        "max_stem_length": max_stem, "max_loop_length": max_loop,
        "dot_bracket": structure,
        "is_stem_loop_like": int(is_stem_loop),
    }


def _count_stems_and_hairpins(structure):
    stems = 0
    in_stem = False
    for c in structure:
        if c == "(":
            if not in_stem:
                stems += 1
                in_stem = True
        else:
            in_stem = False
    hairpin_count = 0
    stack = []
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            if stack:
                open_pos = stack.pop()
                between = structure[open_pos + 1:i]
                if all(x == "." for x in between) and len(between) > 0:
                    hairpin_count += 1
    return stems, hairpin_count


def _max_stem_and_loop(structure):
    stem_runs = re.findall(r"[()]+", structure)
    max_stem = max((len(r) for r in stem_runs), default=0)
    loop_runs = re.findall(r"\.+", structure)
    max_loop = max((len(r) for r in loop_runs), default=0)
    return max_stem, max_loop


def _simplified_structure_features(seq, na_type="DNA"):
    length = len(seq)
    complement = _COMPLEMENT_RNA if na_type.upper() == "RNA" else _COMPLEMENT_DNA
    structure = list("." * length)
    min_loop = 4
    max_len_for_dp = 200

    if length <= max_len_for_dp:
        n_paired, pairs = _nussinov(seq, complement, min_loop)
        for i, j in pairs:
            structure[i] = "("
            structure[j] = ")"
    else:
        n_paired = _greedy_pair_count(seq, complement, min_loop)
        pairs = []

    structure_str = "".join(structure)
    n_paired_total = structure_str.count("(") + structure_str.count(")")
    paired_fraction = n_paired_total / length if length > 0 else 0.0

    mfe_estimate = 0.0
    for i, j in pairs:
        b1, b2 = seq[i], seq[j]
        if b1 in "GC" and b2 in "GC":
            mfe_estimate -= 3.0
        else:
            mfe_estimate -= 2.0

    stems, hairpins = _count_stems_and_hairpins(structure_str)
    max_stem, max_loop = _max_stem_and_loop(structure_str)

    is_stem_loop = (
        paired_fraction >= STEM_LOOP_PAIRED_FRACTION_THRESHOLD
        and hairpins >= STEM_LOOP_MIN_HAIRPIN
    )
    return {
        "mfe": round(mfe_estimate, 2),
        "paired_base_fraction": round(paired_fraction, 4),
        "n_stems": stems, "n_hairpin_loops": hairpins,
        "max_stem_length": max_stem, "max_loop_length": max_loop,
        "dot_bracket": structure_str,
        "is_stem_loop_like": int(is_stem_loop),
    }


def _nussinov(seq, complement, min_loop=4):
    n = len(seq)
    dp = [[0] * n for _ in range(n)]
    for span in range(min_loop + 1, n):
        for i in range(n - span):
            j = i + span
            dp[i][j] = dp[i + 1][j]
            dp[i][j] = max(dp[i][j], dp[i][j - 1])
            if _can_pair(seq[i], seq[j], complement):
                val = (dp[i + 1][j - 1] if i + 1 <= j - 1 else 0) + 1
                dp[i][j] = max(dp[i][j], val)
            for k in range(i + 1, j):
                dp[i][j] = max(dp[i][j], dp[i][k] + dp[k + 1][j])
    pairs = []
    _traceback(dp, seq, complement, 0, n - 1, min_loop, pairs)
    return dp[0][n - 1], pairs


def _traceback(dp, seq, complement, i, j, min_loop, pairs):
    if i >= j:
        return
    if dp[i][j] == dp[i + 1][j]:
        _traceback(dp, seq, complement, i + 1, j, min_loop, pairs)
    elif dp[i][j] == dp[i][j - 1]:
        _traceback(dp, seq, complement, i, j - 1, min_loop, pairs)
    elif _can_pair(seq[i], seq[j], complement) and j - i > min_loop:
        inner = dp[i + 1][j - 1] if i + 1 <= j - 1 else 0
        if dp[i][j] == inner + 1:
            pairs.append((i, j))
            _traceback(dp, seq, complement, i + 1, j - 1, min_loop, pairs)
            return
    for k in range(i + 1, j):
        if dp[i][j] == dp[i][k] + dp[k + 1][j]:
            _traceback(dp, seq, complement, i, k, min_loop, pairs)
            _traceback(dp, seq, complement, k + 1, j, min_loop, pairs)
            return


def _can_pair(b1, b2, complement):
    return complement.get(b1) == b2 or (b1 == "G" and b2 == "U") or (b1 == "U" and b2 == "G")


def _greedy_pair_count(seq, complement, min_loop=4):
    n = len(seq)
    paired = set()
    count = 0
    for gap in range(min_loop + 1, n):
        for i in range(n - gap):
            j = i + gap
            if i not in paired and j not in paired:
                if _can_pair(seq[i], seq[j], complement):
                    paired.add(i)
                    paired.add(j)
                    count += 1
    return count


def _empty_stem_loop_features():
    return {
        "mfe": 0.0, "paired_base_fraction": 0.0, "n_stems": 0,
        "n_hairpin_loops": 0, "max_stem_length": 0, "max_loop_length": 0,
        "dot_bracket": "", "is_stem_loop_like": 0,
    }


# ================================================================
# 最終構造グループ分類
# ================================================================

def classify_structure_group(is_g4_like, is_stem_loop_like):
    if is_g4_like and is_stem_loop_like:
        return "G4_and_stem_loop_like"
    elif is_g4_like:
        return "G4_like"
    elif is_stem_loop_like:
        return "stem_loop_like"
    return "other_or_uncertain"


# ================================================================
# メイン計算関数
# ================================================================

def compute_all_features(df):
    if "sequence" not in df.columns:
        warnings.warn("'sequence' カラムがありません。")
        return df
    logger.info(f"アプタマー特徴量計算開始: {len(df)} 配列")
    all_features = []
    for idx, row in df.iterrows():
        seq = str(row.get("sequence", "")).strip()
        na_type = str(row.get("nucleic_acid_type", "")).strip()
        if not seq:
            features = {}
            features.update(_empty_basic_features())
            features.update(_empty_g4_features())
            features.update(_empty_stem_loop_features())
            features["structure_group"] = "other_or_uncertain"
            all_features.append(features)
            continue
        try:
            basic = compute_basic_features(seq, na_type)
            kmer = compute_kmer_features(seq, k_values=(1, 2, 3))
            g4 = compute_g4_features(seq)
            detected_na = basic.get("nucleic_acid_type_detected", na_type or "DNA")
            stem = compute_stem_loop_features(seq, detected_na)
            structure_group = classify_structure_group(
                bool(g4.get("is_g4_like", 0)),
                bool(stem.get("is_stem_loop_like", 0)),
            )
            features = {}
            features.update(basic)
            features.update(kmer)
            features.update(g4)
            features.update(stem)
            features["structure_group"] = structure_group
        except Exception as e:
            warnings.warn(f"行 {idx} の特徴量計算エラー: {e}")
            features = {}
            features.update(_empty_basic_features())
            features.update(_empty_g4_features())
            features.update(_empty_stem_loop_features())
            features["structure_group"] = "other_or_uncertain"
        all_features.append(features)

    features_df = pd.DataFrame(all_features)
    id_cols = [c for c in ["aptamer_id", "aptamer_name", "sequence", "target_name",
                           "target_name_normalized", "nucleic_acid_type"] if c in df.columns]
    result = pd.concat([df[id_cols].reset_index(drop=True), features_df], axis=1)

    group_counts = result["structure_group"].value_counts()
    logger.info("==== 構造グループ分布 ====")
    for grp, cnt in group_counts.items():
        logger.info(f"  {grp}: {cnt}")
    return result


def main():
    parser = argparse.ArgumentParser(description="アプタマー配列特徴量計算")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    df = safe_read_csv(args.input)
    result = compute_all_features(df)
    safe_write_csv(result, args.output)
    logger.info(f"出力完了: {args.output}  ({len(result)} 行)")


if __name__ == "__main__":
    main()
