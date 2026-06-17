from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_validate, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from src.models.evaluate import classification_metrics, feature_importance_table
from src.utils.io import ensure_parent, read_csv, write_csv
from src.utils.logging import get_logger

LOGGER = get_logger(__name__)
MIN_ROWS = 30
MIN_CLASSES = 2
MIN_PER_CLASS = 5


def feature_columns(df: pd.DataFrame) -> list[str]:
    excluded = {
        "aptamer_id", "aptamer_name", "sequence", "nucleic_acid_type", "target_name", "target_type",
        "kd_unit", "source_database", "reference_doi", "predicted_dot_bracket", "protein_sequence",
        "annotation_text", "uniprot_accession", "uniprot_protein_name", "g4_group", "stem_loop_group",
        "aptamer_structure_group", "charge_group", "hydrophobicity_group", "pI_group", "functional_group",
    }
    numeric = df.select_dtypes(include=[np.number, "bool"]).columns
    return [c for c in numeric if c not in excluded]


def can_train(df: pd.DataFrame, target: str) -> tuple[bool, str]:
    if target not in df.columns:
        return False, f"target column {target} not found"
    usable = df[target].dropna()
    usable = usable[usable != "unknown"]
    counts = usable.value_counts()
    if len(usable) < MIN_ROWS:
        return False, f"need at least {MIN_ROWS} labeled rows, got {len(usable)}"
    if counts.size < MIN_CLASSES:
        return False, "need at least two target classes"
    if counts.min() < MIN_PER_CLASS:
        return False, f"need at least {MIN_PER_CLASS} rows per class, got {counts.to_dict()}"
    return True, "ok"


def train_models(df: pd.DataFrame, target: str) -> tuple[Pipeline, dict[str, object], pd.DataFrame]:
    cols = feature_columns(df)
    data = df[df[target].notna() & (df[target] != "unknown")].copy()
    X = data[cols]
    y = data[target].astype(str)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42, stratify=y)

    candidates = {
        "logistic_regression": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(max_iter=2000, class_weight="balanced")),
        ]),
        "random_forest": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", RandomForestClassifier(n_estimators=300, random_state=42, class_weight="balanced")),
        ]),
        "gradient_boosting": Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("model", GradientBoostingClassifier(random_state=42)),
        ]),
    }
    cv = StratifiedKFold(n_splits=min(5, y.value_counts().min()), shuffle=True, random_state=42)
    summaries = []
    best_name = None
    best_score = -np.inf
    best_model = None
    for name, model in candidates.items():
        scores = cross_validate(
            model,
            X_train,
            y_train,
            cv=cv,
            scoring=["accuracy", "balanced_accuracy", "f1_macro"],
            error_score="raise",
        )
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        metrics = classification_metrics(y_test, pred)
        summary = {
            "model": name,
            "cv_accuracy_mean": float(scores["test_accuracy"].mean()),
            "cv_balanced_accuracy_mean": float(scores["test_balanced_accuracy"].mean()),
            "cv_macro_f1_mean": float(scores["test_f1_macro"].mean()),
            **metrics,
        }
        summaries.append(summary)
        if summary["cv_macro_f1_mean"] > best_score:
            best_score = summary["cv_macro_f1_mean"]
            best_name = name
            best_model = model
    assert best_model is not None
    best_model.fit(X, y)
    model_step = best_model.named_steps["model"]
    if hasattr(model_step, "feature_importances_"):
        importance = feature_importance_table(cols, model_step.feature_importances_)
    elif hasattr(model_step, "coef_"):
        importance = feature_importance_table(cols, np.mean(np.abs(model_step.coef_), axis=0))
    else:
        importance = pd.DataFrame({"feature": cols, "importance": np.nan})
    metadata = {"target": target, "features": cols, "best_model": best_name, "summaries": summaries}
    return best_model, metadata, importance


def markdown_summary(metadata: dict[str, object]) -> str:
    summaries = pd.DataFrame(metadata["summaries"])
    keep = ["model", "cv_accuracy_mean", "cv_balanced_accuracy_mean", "cv_macro_f1_mean", "accuracy", "balanced_accuracy", "macro_f1"]
    table = dataframe_to_markdown(summaries[keep])
    return (
        f"目的変数 `{metadata['target']}` に対して `{metadata['best_model']}` を最良モデルとして保存した。\n\n"
        + table
        + "\n\nこの結果は探索的解析であり、実験的検証が必要である。"
    )


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return "該当なし"
    headers = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in df.columns) + " |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/processed/merged_features.csv")
    parser.add_argument("--target", default="hydrophobicity_group")
    parser.add_argument("--model-output", default="models/baseline_model.pkl")
    parser.add_argument("--metrics-output", default="data/processed/baseline_model_metrics.csv")
    parser.add_argument("--importance-output", default="data/processed/baseline_feature_importance.csv")
    parser.add_argument("--summary-output", default="data/processed/baseline_model_summary.md")
    args = parser.parse_args()
    df = read_csv(args.input)
    ok, reason = can_train(df, args.target)
    if not ok:
        msg = f"Baseline model skipped: {reason}. This positive-only dataset should not be overinterpreted."
        LOGGER.warning(msg)
        Path(ensure_parent(args.summary_output)).write_text(msg, encoding="utf-8")
        return
    model, metadata, importance = train_models(df, args.target)
    joblib.dump({"pipeline": model, "metadata": metadata}, ensure_parent(args.model_output))
    write_csv(pd.DataFrame(metadata["summaries"]), args.metrics_output)
    write_csv(importance, args.importance_output)
    Path(ensure_parent(args.summary_output)).write_text(markdown_summary(metadata), encoding="utf-8")
    LOGGER.info("Saved baseline model to %s", args.model_output)


if __name__ == "__main__":
    main()
