#!/usr/bin/env python3
"""Extract fixed-length PAM counts from gzipped FASTQ files."""

from __future__ import annotations

import argparse
import csv
import gzip
import pickle
from pathlib import Path

from Bio import SeqIO


VALID_BASES = set("ACGT")


def extract_counts(fastq: Path, five_seed: str, three_seed: str, pam_len: int) -> tuple[dict[str, int], int]:
    counts: dict[str, int] = {}
    total = 0
    with gzip.open(fastq, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)
            five_idx = seq.find(five_seed)
            three_idx = seq.find(three_seed)
            if five_idx == -1 or three_idx == -1:
                continue
            pam = seq[five_idx + len(five_seed) : three_idx]
            if len(pam) != pam_len:
                continue
            if any(base not in VALID_BASES for base in pam):
                continue
            counts[pam] = counts.get(pam, 0) + 1
            total += 1
    return counts, total


def read_sample_sheet(path: Path) -> list[tuple[str, Path]]:
    rows: list[tuple[str, Path]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"sample", "fastq"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"sample sheet missing columns: {sorted(missing)}")
        for row in reader:
            rows.append((row["sample"], Path(row["fastq"])))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--five-seed", default="CTAGAGGATC")
    parser.add_argument("--three-seed", default="AGGGCGACAC")
    parser.add_argument("--pam-len", default=5, type=int)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for sample, fastq in read_sample_sheet(args.sample_sheet):
        if not fastq.exists():
            raise FileNotFoundError(fastq)
        counts, total = extract_counts(fastq, args.five_seed, args.three_seed, args.pam_len)
        out_path = args.out_dir / f"{sample}_counts.p"
        with out_path.open("wb") as handle:
            pickle.dump((counts, total), handle)
        covered = sum(1 for value in counts.values() if value > 0)
        print(f"{sample}\tvalid_reads={total}\tcovered_pams={covered}\toutput={out_path}")


if __name__ == "__main__":
    main()
