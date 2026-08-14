#!/usr/bin/env python3
"""
plot_sasa.py

Plots SASA vs time and a probability-density distribution for one or more
replicas, following the same calling convention as plot_rg.py:

    python3 plot_sasa.py --outdir <dir> --label <rep_label> file1.xvg [file2.xvg ...]

Each input .xvg is expected to have whitespace-separated numeric columns:
time in column 1, SASA area in --column (1-indexed, default 2 = gmx sasa's
default total-area column).

Distribution bin width is auto-selected via the Freedman-Diaconis rule
(2 * IQR / n^(1/3)) — robust to outliers/skew and standard for physical
time-series data — unless --bin-width is given explicitly to override it.
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless — these scripts run over SSH on remote machines
import matplotlib.pyplot as plt


def read_xvg(path, column):
    time, values = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) < column:
                continue
            try:
                t = float(parts[0])
                v = float(parts[column - 1])
            except ValueError:
                continue
            time.append(t)
            values.append(v)
    return np.array(time), np.array(values)


def freedman_diaconis_bin_width(data):
    q75, q25 = np.percentile(data, [75, 25])
    iqr = q75 - q25
    n = len(data)
    if iqr == 0 or n == 0:
        return None  # degenerate data — fall back to numpy's own default rule
    return 2 * iqr / (n ** (1 / 3))


def main():
    ap = argparse.ArgumentParser(description="Plot SASA vs time + distribution")
    ap.add_argument("xvg_files", nargs="+", help="gmx sasa .xvg output file(s)")
    ap.add_argument("--outdir", default=".", help="output directory")
    ap.add_argument("--label", required=True, help="replica label, e.g. rep1")
    ap.add_argument("--column", type=int, default=2,
                     help="1-indexed column to read as SASA (default: 2, gmx sasa's total-area column)")
    ap.add_argument("--bin-width", type=float, default=None,
                     help="override the automatic Freedman-Diaconis bin width")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for xvg in args.xvg_files:
        base = os.path.splitext(os.path.basename(xvg))[0]
        time, sasa = read_xvg(xvg, args.column)
        if len(sasa) == 0:
            print(f"WARNING: no numeric data parsed from {xvg} — skipping.")
            continue

        # ── Time series plot ──────────────────────────────────────────────
        plt.figure(figsize=(8, 6))
        plt.plot(time, sasa, linewidth=0.8)
        plt.xlabel("Time (ps)")
        plt.ylabel("SASA (nm²)")
        plt.title(f"SASA vs time — {args.label}")
        plt.tight_layout()
        ts_png = os.path.join(args.outdir, f"{base}_vs_time.png")
        plt.savefig(ts_png, dpi=300)
        plt.close()

        # ── Stats ──────────────────────────────────────────────────────────
        avg = np.mean(sasa)
        std = np.std(sasa)
        stats_file = os.path.join(args.outdir, f"{base}_average_stdev.dat")
        with open(stats_file, "w") as f:
            f.write("# Statistics of SASA\n")
            f.write(f"Number of points: {len(sasa)}\n")
            f.write(f"Average: {avg:.6f}\n")
            f.write(f"Standard Deviation: {std:.6f}\n")

        # ── Distribution (auto bin width via Freedman-Diaconis) ────────────
        bin_width = args.bin_width if args.bin_width else freedman_diaconis_bin_width(sasa)
        if bin_width is None or bin_width <= 0:
            bin_edges_arg = "auto"  # numpy fallback for degenerate/tiny samples
        else:
            bin_edges_arg = np.arange(sasa.min(), sasa.max() + bin_width, bin_width)

        hist, edges = np.histogram(sasa, bins=bin_edges_arg, density=True)
        bin_centers = 0.5 * (edges[:-1] + edges[1:])
        actual_width = edges[1] - edges[0]

        dist_dat = os.path.join(args.outdir, f"{base}_Pdist.dat")
        with open(dist_dat, "w") as f:
            source = "manual override" if args.bin_width else "auto (Freedman-Diaconis)"
            f.write(f"# Bin_center    Probability_density   (bin width = {actual_width:.6f}, {source})\n")
            for center, prob in zip(bin_centers, hist):
                f.write(f"{center:.6f}    {prob:.6f}\n")

        plt.figure(figsize=(7, 4))
        plt.plot(bin_centers, hist, drawstyle="steps-mid")
        plt.xlabel("SASA (nm²)")
        plt.ylabel("Probability Density")
        plt.title(f"SASA distribution — {args.label} (bin width = {actual_width:.3f})")
        plt.tight_layout()
        dist_png = os.path.join(args.outdir, f"{base}_Pdist.png")
        plt.savefig(dist_png, dpi=300)
        plt.close()

        print(f"→ {base}: avg={avg:.4f} nm^2, std={std:.4f} nm^2, bin_width={actual_width:.4f}")
        print(f"   time series : {ts_png}")
        print(f"   distribution: {dist_png}  ({dist_dat})")


if __name__ == "__main__":
    main()