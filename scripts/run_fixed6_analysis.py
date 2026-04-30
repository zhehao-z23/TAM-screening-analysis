#!/usr/bin/env python3
"""Run the fixed six-configuration TAM/PAM screen analysis."""

from __future__ import annotations

import argparse
import itertools
import pickle
from pathlib import Path

import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def pam_keys(pam_len: int) -> list[str]:
    return ["".join(p) for p in itertools.product("ATGC", repeat=pam_len)]


def load_pct_vector(counts_dir: Path, label: str, keys: list[str]) -> tuple[int, int, np.ndarray]:
    with (counts_dir / f"{label}_counts.p").open("rb") as handle:
        counts, total = pickle.load(handle)
    vec = np.array([(counts.get(k, 0) / total * 100) if total > 0 else 0 for k in keys], dtype=float)
    covered = sum(1 for key in keys if counts.get(key, 0) > 0)
    return total, covered, vec


def bits_matrix(seqs: list[str], pam_len: int) -> pd.DataFrame | None:
    if not seqs:
        return None
    counts = pd.DataFrame(0, index=range(pam_len), columns=["A", "C", "G", "T"])
    for seq in seqs:
        for idx, base in enumerate(seq):
            if base in counts.columns:
                counts.loc[idx, base] += 1
    probs = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    entropy = -(probs * np.log2(probs + 1e-9)).sum(axis=1)
    return probs.multiply(2 - entropy, axis=0)


def draw_logo(matrix: pd.DataFrame | None, ax: plt.Axes, title: str) -> None:
    if matrix is None:
        ax.text(0.5, 0.5, "No motif", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return
    logomaker.Logo(matrix, color_scheme="classic", ax=ax)
    ax.set_ylim(0, 2)
    ax.set_title(title)


def lower_threshold(fc: np.ndarray, rule: str) -> tuple[float, float, float]:
    log_fc = np.log(fc + 1e-30)
    mu = float(np.mean(log_fc))
    sd = float(np.std(log_fc))
    if rule == "sigma1.5":
        lower = float(np.exp(mu - 1.5 * sd))
    elif rule == "sigma5":
        lower = float(np.exp(mu - 5.0 * sd))
    elif rule == "fc0.15":
        lower = 0.15
    else:
        raise ValueError(f"unknown rule: {rule}")
    return lower, mu, sd


def upper_threshold(rule: str, mu: float, sd: float) -> float:
    if rule == "sigma1.5":
        return float(np.exp(mu + 1.5 * sd))
    if rule == "sigma5":
        return float(np.exp(mu + 5.0 * sd))
    if rule == "fc0.15":
        return 1 / 0.15
    raise ValueError(f"unknown rule: {rule}")


def run_sample(args: argparse.Namespace, keys: list[str], sample: str, x_pct: np.ndarray) -> list[dict[str, object]]:
    total, covered, y_pct = load_pct_vector(args.counts_dir, sample, keys)
    print(f"[sample] {sample}: total_reads={total}, covered={covered}/{len(keys)}")

    raw_fc = y_pct / (x_pct + args.epsilon)
    correction_factor = float(np.median(raw_fc))
    y_corr = y_pct / correction_factor
    corr_fc = y_corr / (x_pct + args.epsilon)

    modes = {
        "raw": {"y": y_pct, "fc": raw_fc, "correction_factor": 1.0},
        "corrected": {"y": y_corr, "fc": corr_fc, "correction_factor": correction_factor},
    }
    rows: list[dict[str, object]] = []
    out_dir = args.out_dir / sample
    out_dir.mkdir(parents=True, exist_ok=True)

    for mode_name, mode in modes.items():
        y_use = mode["y"]
        fc_use = mode["fc"]
        for rule in ["sigma1.5", "sigma5", "fc0.15"]:
            lower, mu, sd = lower_threshold(fc_use, rule)
            upper = upper_threshold(rule, mu, sd)
            depleted = fc_use <= lower
            enriched = fc_use >= upper

            dep_idx = np.where(depleted)[0]
            enr_idx = np.where(enriched)[0]
            dep_top = dep_idx[np.argsort(fc_use[dep_idx])][: args.top_n] if len(dep_idx) else []
            enr_top = enr_idx[np.argsort(fc_use[enr_idx])[::-1]][: args.top_n] if len(enr_idx) else []

            fig, axs = plt.subplots(2, 4, figsize=(22, 10))
            x_line = np.logspace(np.log10(args.global_min), np.log10(args.global_max), 200)
            log2fc = np.log2(fc_use + 1e-30)

            panels = [
                (0, depleted, lower, "red", "Depleted", dep_idx, dep_top),
                (1, enriched, upper, "blue", "Enriched", enr_idx, enr_top),
            ]
            for row, mask, cutoff, color, label, idx, top_idx in panels:
                ax = axs[row, 0]
                ax.scatter(x_pct, y_use, s=10, color="#d3d3d3", alpha=0.35, edgecolor="none")
                ax.scatter(x_pct[mask], y_use[mask], s=24, color=color, alpha=0.85, edgecolor="none", label=label)
                ax.plot(x_line, x_line, "k-", lw=1, alpha=0.6)
                ax.plot(x_line, x_line * cutoff, "--", color=color, lw=1, alpha=0.7, label=f"cutoff={cutoff:.4g}")
                ax.set_xscale("log"); ax.set_yscale("log")
                ax.set_xlim(args.global_min, args.global_max); ax.set_ylim(args.global_min, args.global_max)
                ax.set_xlabel("Control abundance (%)")
                ax.set_ylabel("Experiment abundance (%)")
                ax.set_title(f"{sample} | {mode_name} | {rule} | {label}")
                ax.legend(frameon=False, loc="upper left")

                ax = axs[row, 1]
                ax.hist(log2fc, bins=40, color="#7f8c8d", alpha=0.85)
                ax.axvline(np.log2(cutoff), color=color, ls="--", lw=1.5)
                ax.set_xlabel("log2(FC)")
                ax.set_ylabel("Count of PAMs")
                ax.set_title(f"FC distribution ({label.lower()} cutoff)")

                seqs = [keys[i] for i in idx]
                top_seqs = [keys[i] for i in top_idx]
                draw_logo(bits_matrix(seqs, args.pam_len), axs[row, 2], f"{label} motif (n={len(seqs)})")
                draw_logo(bits_matrix(top_seqs, args.pam_len), axs[row, 3], f"Top{min(args.top_n, len(top_seqs))} {label.lower()} motif")

            fig.tight_layout()
            plot_file = out_dir / f"{sample}_{mode_name}_{rule}.png"
            fig.savefig(plot_file, dpi=180)
            plt.close(fig)

            rows.append({
                "experiment": sample,
                "mode": mode_name,
                "rule": rule,
                "correction_factor": mode["correction_factor"],
                "threshold_lower": lower,
                "threshold_upper": upper,
                "mu_logfc": mu,
                "sd_logfc": sd,
                "depleted_count": int(np.sum(depleted)),
                "enriched_count": int(np.sum(enriched)),
                "total_pams": len(keys),
                "depleted_fraction": float(np.sum(depleted) / len(keys)),
                "enriched_fraction": float(np.sum(enriched) / len(keys)),
                "plot_file": str(plot_file),
            })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-dir", required=True, type=Path)
    parser.add_argument("--control", required=True)
    parser.add_argument("--samples", required=True, nargs="+")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--pam-len", default=5, type=int)
    parser.add_argument("--epsilon", default=1e-9, type=float)
    parser.add_argument("--global-min", default=1e-3, type=float)
    parser.add_argument("--global-max", default=2.0, type=float)
    parser.add_argument("--top-n", default=20, type=int)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    keys = pam_keys(args.pam_len)
    total, covered, x_pct = load_pct_vector(args.counts_dir, args.control, keys)
    print(f"[control] {args.control}: total_reads={total}, covered={covered}/{len(keys)}")

    rows: list[dict[str, object]] = []
    for sample in args.samples:
        rows.extend(run_sample(args, keys, sample, x_pct))

    summary = pd.DataFrame(rows)
    summary_file = args.out_dir / "summary_all_split_6configs.csv"
    summary.to_csv(summary_file, index=False)
    print(f"wrote {summary_file}")


if __name__ == "__main__":
    main()
