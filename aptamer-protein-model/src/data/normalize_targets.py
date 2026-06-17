from __future__ import annotations

import re


def normalize_target_name(name: object) -> str:
    if name is None:
        return ""
    value = str(name).strip()
    value = re.sub(r"\s+", " ", value)
    value = re.sub(r"\s*\([^)]*\)\s*$", "", value)
    return value


def normalize_target_type(value: object) -> str:
    text = "" if value is None else str(value).strip().lower()
    if text in {"protein", "proteins", "peptide"}:
        return "protein"
    return text or "unknown"
