#!/usr/bin/env python3
"""
plot_dssp.py

Parses a DSSP .dat file (auto-detected in the current directory if not specified)
and computes per-residue secondary structure distributions for alpha-helix,
beta-sheet, turn, and coil.

Usage:
    python3 plot_dssp.py
    python3 plot_dssp.py DSSP_Protein_rep1.dat
"""

import argparse
import glob
import os
import re
import warnings

# Suppress harmless Matplotlib UserWarnings (e.g. Axes3D import warnings)
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

import numpy as np
import matplotlib
matplotlib.use("Agg")  # Headless mode for remote clusters/SSH
import matplotlib.pyplot as plt

# Publication-ready figure styling
plt.rcParams.update({
    "font.family": "serif",
    "axes.linewidth": 2,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 2,
    "ytick.major.width": 2,
    "legend.frameon": False,
})


def find_target_dat_file(outdir="."):
    """Auto-detect the generated DSSP .dat file in the target directory."""
    # Priority 1: DSSP_Protein_*.dat or DSSP_*.dat files
    candidates = glob.glob(os.path.join(outdir, "DSSP_*.dat"))
    
    # Priority 2: Any .dat file excluding previously generated output distributions
    if not candidates:
        all_dat = glob.glob(os.path.join(outdir, "*.dat"))
        ignore_list = ["_distribution_", "_summary.dat"]
        candidates = [f for f in all_dat if not any(ign in f for ign in ignore_list)]

    if not candidates:
        raise FileNotFoundError(
            f"No valid DSSP .dat input file found in '{outdir}'. "
            "Please run script_dssp_replica.sh first or pass the file path directly."
        )

    # Return the first match (or DSSP_Protein_*.dat if present)
    candidates.sort()
    return candidates[0]


def infer_replica_label(filepath):
    """Infer replica label from filename (e.g. DSSP_Protein_rep1.dat -> rep1)."""
    basename = os.path.splitext(os.path.basename(filepath))[0]
    
    # Match rep1, R1, replica1, etc. from filename
    match = re.search(r'(rep\d+|R\d+|replica\d+)', basename, re.IGNORECASE)
    if match:
        return match.group(1)
    
    # Fallback: Extract whatever comes after the last underscore
    parts = basename.split("_")
    if len(parts) > 1:
        return parts[-1]
    
    # Ultimate fallback: Current folder name
    return os.path.basename(os.getcwd())


def main():
    ap = argparse.ArgumentParser(description="Per-residue targeted DSSP distribution from .dat")
    ap.add_argument("dat_file", nargs="?", default=None, help="Input .dat file (optional, auto-detected if omitted)")
    ap.add_argument("--outdir", default=".", help="Output directory (default: current directory)")
    ap.add_argument("--label", default=None, help="Replica label (optional, auto-inferred if omitted)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Resolve input file
    input_file = args.dat_file if args.dat_file else find_target_dat_file(args.outdir)
    
    # 2. Resolve label
    label = args.label if args.label else infer_replica_label(input_file)

    print(f"--> Input file    : {input_file}")
    print(f"--> Replica Label : {label}")

    # ── Parse the .dat file ───────────────────────────────────────────────
    with open(input_file, "r") as f:
        lines = f.read().splitlines()

    data_lines = []
    for line in lines:
        s = line.strip()
        # Skip empty lines or header integers (like residue count line '42')
        if not s or s.isdigit():
            continue
        data_lines.append(list(s))

    if not data_lines:
        raise ValueError(f"No valid frame sequences found in {input_file}")

    # Matrix shape: (n_residues, n_frames)
    matrix = np.array(data_lines).T 
    n_res, n_frames = matrix.shape
    residues = np.arange(1, n_res + 1)

    print(f"--> Dataset stats : {n_res} residues x {n_frames} frames")

    # ── Target SS Types ───────────────────────────────────────────────────
    SS_TARGETS = {
        "alpha_helix": {
            "chars": ["H"], 
            "file": "H_distribution_alpha_helix.dat", 
            "label": r"$\alpha$-Helix", 
            "color": "#D55E00"
        },
        "beta_sheet": {
            "chars": ["E"], 
            "file": "E_distribution_beta_sheet.dat", 
            "label": r"$\beta$-Sheet", 
            "color": "#0072B2"
        },
        "turn": {
            "chars": ["T"], 
            "file": "T_distribution_turn.dat", 
            "label": "Turn", 
            "color": "#a056bf"
        },
        "coil": {
            "chars": ["~", "C", "-"], 
            "file": "coil_distribution_coil.dat", 
            "label": "Coil", 
            "color": "#999999"
        },
    }

    # ── Compute & Save Distributions ──────────────────────────────────────
    dist_data = {}
    
    for target_key, info in SS_TARGETS.items():
        frac = np.zeros(n_res)
        for ch in info["chars"]:
            frac += (matrix == ch).sum(axis=1) / n_frames
        
        dist_data[target_key] = frac
        
        # Save per-residue profile directly using the clean standard filename
        out_file = os.path.join(args.outdir, info['file'])
        with open(out_file, "w") as f:
            for i, res in enumerate(residues):
                f.write(f"{res:4d}    {frac[i]:.6f}\n")
        
        print(f"    Saved: {out_file}")

    # ── Plot Distributions ────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(8, n_res * 0.15), 5))
    fig.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15)
    
    for target_key, info in SS_TARGETS.items():
        ax.plot(residues, dist_data[target_key], linewidth=2.5,
                color=info["color"], label=info["label"])
        ax.fill_between(residues, dist_data[target_key], color=info["color"], alpha=0.1)

    ax.set_xlabel("Residue", fontsize=14)
    ax.set_ylabel("SS Probability", fontsize=14)
    ax.set_title(f"Secondary Structure Distributions — {label}", fontsize=16)
    ax.legend(loc="upper right", fontsize=12)
    
    ax.set_xlim(residues[0], residues[-1])
    ax.set_ylim(0, 1)
    ax.tick_params(axis="both", width=2, length=8, labelsize=12)
    
    dist_png = os.path.join(args.outdir, f"{label}_ss_distributions.png")
    plt.savefig(dist_png, dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"--> Plot saved  : {dist_png}\nDone.")

if __name__ == "__main__":
    main()