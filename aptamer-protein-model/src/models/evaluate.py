from __future__ import annotations

import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, f1_score


def classification_metrics(y_true, y_pred) -> dict[str, object]:
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "macro_f1": f1_score(y_true, y_pred, average="macro", zero_division=0),
        "confusion_matrix": confusion_matrix(y_true, y_pred).tolist(),
    }


def feature_importance_table(feature_names: list[str], importances) -> pd.DataFrame:
    return pd.DataFrame({"feature": feature_names, "importance": importances}).sort_values("importance", ascending=False)
