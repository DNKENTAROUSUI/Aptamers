from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def save_bar_counts(series: pd.Series, title: str, output: str | Path) -> None:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    counts = series.fillna("unknown").value_counts()
    fig, ax = plt.subplots(figsize=(8, 4.8))
    counts.plot(kind="bar", ax=ax, color="#4C78A8")
    ax.set_title(title)
    ax.set_ylabel("count")
    ax.set_xlabel("")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def save_hist(series: pd.Series, title: str, output: str | Path, bins: int = 30) -> None:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    values = pd.to_numeric(series, errors="coerce").dropna()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if values.empty:
        ax.text(0.5, 0.5, "No numeric data", ha="center", va="center")
    else:
        ax.hist(values, bins=min(bins, max(5, len(values))), color="#59A14F", edgecolor="white")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def save_heatmap(table: pd.DataFrame, title: str, output: str | Path) -> None:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    fig_width = max(9, 1.2 * max(1, len(table.columns)) + 4)
    fig_height = max(5.5, 0.65 * max(1, len(table.index)) + 2.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    if table.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
    else:
        values = table.values.astype(float)
        image = ax.imshow(values, cmap="Blues", aspect="auto", vmin=0, vmax=max(1.0, float(np.nanmax(values))))
        ax.set_xticks(range(len(table.columns)))
        ax.set_xticklabels(table.columns, rotation=35, ha="right")
        ax.set_yticks(range(len(table.index)))
        ax.set_yticklabels(table.index)
        threshold = float(np.nanmax(values)) * 0.55 if values.size else 0.0
        for i in range(table.shape[0]):
            for j in range(table.shape[1]):
                value = float(table.iloc[i, j])
                color = "white" if value > threshold else "#1f2933"
                ax.text(j, i, int(value), ha="center", va="center", color=color, fontweight="bold")
        fig.colorbar(image, ax=ax, label="count")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)
