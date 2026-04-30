#!/usr/bin/env python3
"""Compare N5 libraries and visualize library bias."""

from __future__ import annotations

import argparse
import itertools
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def pam_keys(pam_len: int) -> list[str]:
    return ["".join(p) for p in itertools.product("ATGC", repeat=pam_len)]


def load_vec(path: Path, keys: list[str]) -> tuple[int, int, np.ndarray]:
    with path.open("rb") as handle:
        counts, total = pickle.load(handle)
    vec = np.array([counts.get(k, 0) / total if total > 0 else 0 for k in keys], dtype=float)
    covered = sum(1 for key in keys if counts.get(key, 0) > 0)
    return total, covered, vec


def gini(x: np.ndarray) -> float:
    values = np.sort(np.asarray(x, dtype=float))
    if values.sum() == 0:
        return 0.0
    n = len(values)
    cum = np.cumsum(values)
    return float((n + 1 - 2 * cum.sum() / cum[-1]) / n)


def jsd(p: np.ndarray, q: np.ndarray, eps: float = 1e-12) -> float:
    p = np.asarray(p, dtype=float) + eps
    q = np.asarray(q, dtype=float) + eps
    p /= p.sum(); q /= q.sum()
    m = (p + q) / 2
    return float(np.sqrt(0.5 * (np.sum(p * np.log2(p / m)) + np.sum(q * np.log2(q / m)))))


def parse_n5(values: list[str]) -> dict[str, Path]:
    parsed: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--n5 entries must be NAME=/path/to/N5_counts.p")
        name, path = value.split("=", 1)
        parsed[name] = Path(path)
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n5", action="append", required=True, help="NAME=/path/to/N5_counts.p")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--pam-len", default=5, type=int)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    keys = pam_keys(args.pam_len)
    paths = parse_n5(args.n5)

    vecs: dict[str, np.ndarray] = {}
    metrics: list[dict[str, object]] = []
    for name, path in paths.items():
        total, covered, vec = load_vec(path, keys)
        vecs[name] = vec
        metrics.append({
            "n5_name": name,
            "total_reads": total,
            "covered_pams": covered,
            "top1_frac": float(vec.max()),
            "top10_frac": float(np.sort(vec)[-10:].sum()),
            "entropy_bits": float(-(vec * np.log2(vec + 1e-12)).sum()),
            "gini": gini(vec),
        })

    names = list(vecs)
    pearson = np.eye(len(names))
    pairs: list[dict[str, object]] = []
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            pearson[i, j] = np.corrcoef(vecs[a], vecs[b])[0, 1]
        for b in names[i + 1:]:
            log2fc = np.log2((vecs[b] + 1e-12) / (vecs[a] + 1e-12))
            pairs.append({
                "pair": f"{a}__vs__{b}",
                "pearson": float(np.corrcoef(vecs[a], vecs[b])[0, 1]),
                "jsd": jsd(vecs[a], vecs[b]),
                "median_log2fc": float(np.median(log2fc)),
                "mad_log2fc": float(np.median(np.abs(log2fc - np.median(log2fc)))),
                "p95_abs_log2fc": float(np.percentile(np.abs(log2fc), 95)),
            })

    metrics_df = pd.DataFrame(metrics)
    metrics_df["avg_pearson_to_others"] = [float((pearson[i].sum() - 1) / (len(names) - 1)) for i in range(len(names))]
    metrics_df.to_csv(args.out_dir / "n5_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(pairs).to_csv(args.out_dir / "n5_pairwise.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    image = ax.imshow(pearson, vmin=0.65, vmax=1.0, cmap="viridis")
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
    for i in range(len(names)):
        for j in range(len(names)):
            ax.text(j, i, f"{pearson[i, j]:.3f}", ha="center", va="center", color="white", fontsize=8)
    ax.set_title("N5 Pearson Correlation")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(args.out_dir / "n5_pearson_heatmap.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    for name, vec in vecs.items():
        ax.hist(np.log10(vec * 100 + 1e-9), bins=50, histtype="step", linewidth=1.6, label=name)
    ax.set_xlabel("log10(PAM abundance %)")
    ax.set_ylabel("PAM count")
    ax.set_title("N5 abundance distribution")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out_dir / "n5_abundance_distribution.png", dpi=180)
    plt.close(fig)

    print(metrics_df.to_string(index=False))
    print(f"wrote {args.out_dir}")


if __name__ == "__main__":
    main()
