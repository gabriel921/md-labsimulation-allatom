#!/usr/bin/env python3
"""
plot_rg.py -- Rg vs time plot from one or more GROMACS gyrate .xvg files.

Usage:
    python3 plot_rg.py --outdir DIR --label REPLICA_LABEL file1.xvg [file2.xvg ...]

Produces: DIR/Rg_vs_time_REPLICA_LABEL.png
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_xvg(path):
    """Return (time, rg) arrays from a GROMACS xvg file, skipping headers."""
    time, rg = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            time.append(float(parts[0]))
            rg.append(float(parts[1]))  # col 2 = Rg; cols 3-5 are Rg_x/y/z
    return time, rg


def group_name_from_filename(path):
    # Expects Rg_<group>_<replabel>.xvg
    stem = Path(path).stem
    m = re.match(r"Rg_(.+)_[^_]+$", stem)
    return m.group(1) if m else stem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--label", default="replica")
    ap.add_argument("xvg_files", nargs="+")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))

    for f in args.xvg_files:
        try:
            time, rg = read_xvg(f)
        except Exception as e:
            print(f"WARNING: could not read {f}: {e}", file=sys.stderr)
            continue
        if not time:
            print(f"WARNING: no data points in {f}", file=sys.stderr)
            continue
        time_ns = [t / 1000.0 for t in time]  # gmx default time unit is ps
        ax.plot(time_ns, rg, label=group_name_from_filename(f), linewidth=1.2)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Radius of gyration (nm)")
    ax.set_title(f"Rg vs time — {args.label}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()

    outfile = outdir / f"Rg_vs_time_{args.label}.png"
    fig.savefig(outfile, dpi=150)
    print(f"Plot written: {outfile}")


if __name__ == "__main__":
    main()
