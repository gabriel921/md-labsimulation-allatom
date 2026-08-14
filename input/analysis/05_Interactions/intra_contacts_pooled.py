#!/usr/bin/env python3
"""
Pooled intra-chain residue-residue contact analysis across all replicas of
the Aβ42 monomer pipeline.

Structured to closely mirror the original chain-pair contacts script:
for each replica, walk the trajectory and compute the total number of
residue pairs in contact per frame, write a per-replica .dat file, and
plot ONE line per replica on a single shared figure.

The only real change from "chain-pair" to "intra-chain" is which residue
pairs get counted: instead of pairs drawn from two different chains, pairs
are drawn from WITHIN each chain named in CHAIN_GROUPS (residue self-pairs
excluded; sequence-neighbor exclusion is optional via RESIDUE_EXCLUDE).

If more than one chain is given in CHAIN_GROUPS, each chain's intra-contacts
are counted separately and also summed into one total_contacts column per
frame (same pattern as the original script summing multiple chain-pairs
into one total) — this is what feeds the single overlaid plot.

All configuration comes from environment variables (set by the calling
script_intracontacts_pooled.sh launcher).

Required env vars:
    TOPOLOGY    - path to nowater.tpr (MDAnalysis topology, has bonds/masses)
    NDXFILE     - path to index_analysis.ndx
    TRAJFILES   - comma-separated list of trajectory paths, one per replica
    REPNAMES    - comma-separated replica labels, same order as TRAJFILES

Optional env vars (defaults shown):
    CHAIN_GROUPS        = Protein   (comma-separated NDX group names; each
                                      chain's intra-contacts are counted
                                      separately and summed into the total)
    CUTOFF              = 0.35      (Å — MDAnalysis distance units. NOTE:
                                      0.35 Å is sub-atomic; if you meant
                                      0.35 nm, set CUTOFF=3.5)
    RESIDUE_EXCLUDE     = 0         (exclude |i-j| <= this many residues
                                      apart; 0 = no skipping, only i==j
                                      itself is excluded)
    ANALYSIS_INTERVAL_NS= 1         (frame stride target, in ns)
    MA_WINDOW           = 0         (moving-average window; 0 = off)
    OUTDIR              = contacts_intra
"""

import sys
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from itertools import combinations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tqdm import tqdm


def env(name, default=None, cast=str):
    val = os.environ.get(name)
    if val is None or val == "":
        return default
    return cast(val)


# ====================== CONFIG (from environment) ======================
topology = env("TOPOLOGY")
ndxfile = env("NDXFILE")
trajfiles_raw = env("TRAJFILES")
repnames_raw = env("REPNAMES")

if not topology:
    sys.exit("ERROR: TOPOLOGY env var must be set (path to nowater.tpr)")
if not ndxfile:
    sys.exit("ERROR: NDXFILE env var must be set (path to index_analysis.ndx)")
if not trajfiles_raw:
    sys.exit("ERROR: TRAJFILES env var must be set (comma-separated trajectory paths)")

trajfiles = [t.strip() for t in trajfiles_raw.split(",") if t.strip()]
repnames = [r.strip() for r in (repnames_raw.split(",") if repnames_raw else [])] or \
           [f"rep{i+1}" for i in range(len(trajfiles))]

if len(repnames) != len(trajfiles):
    sys.exit(f"ERROR: REPNAMES ({len(repnames)}) and TRAJFILES ({len(trajfiles)}) count mismatch")

chain_groups = [c.strip() for c in env("CHAIN_GROUPS", "Protein").split(",") if c.strip()]
cutoff = env("CUTOFF", 0.35, float)
residue_exclude = env("RESIDUE_EXCLUDE", 0, int)
analysis_interval_ns = env("ANALYSIS_INTERVAL_NS", 0, float)
ma_window = env("MA_WINDOW", 0, int)
outdir = env("OUTDIR", "contacts_intra")

os.makedirs(outdir, exist_ok=True)

print("=== Pooled intra-chain contacts — configuration ===")
print(f"  TOPOLOGY              : {topology}")
print(f"  NDXFILE                : {ndxfile}")
print(f"  Chain groups ({len(chain_groups)})    : {chain_groups}")
print(f"  Replicas ({len(trajfiles)})           :")
for r, t in zip(repnames, trajfiles):
    print(f"    - {r}: {t}")
print(f"  CUTOFF                 : {cutoff} Å"
      + ("  <-- NOTE: sub-atomic if you meant nm; use CUTOFF=3.5 for 0.35 nm" if cutoff < 1.0 else ""))
print(f"  RESIDUE_EXCLUDE        : {residue_exclude}  (|i-j| <= this many residues apart is skipped)")
print(f"  ANALYSIS_INTERVAL_NS   : {analysis_interval_ns} ns")
print(f"  MA_WINDOW              : {ma_window}")
print(f"  OUTDIR                 : {outdir}")


# ====================== NDX PARSING ======================
def read_ndx(ndx_file):
    groups = {}
    name = None
    indices = []
    with open(ndx_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("["):
                if name is not None:
                    groups[name] = indices
                name = line.replace("[", "").replace("]", "").strip()
                indices = []
            elif line and not line.startswith(";"):
                indices.extend([int(x) - 1 for x in line.split()])
    if name is not None:
        groups[name] = indices
    return groups


try:
    ndx = read_ndx(ndxfile)
except IOError:
    raise FileNotFoundError(f"Could not find or open the index file: {ndxfile}")

for cg in chain_groups:
    if cg not in ndx:
        raise KeyError(
            f"Chain group '{cg}' not found in {ndxfile}. "
            f"Available groups: {list(ndx.keys())}"
        )


def moving_average(x, window):
    if window <= 1:
        return x
    pad_left = window // 2
    pad_right = window - 1 - pad_left
    xpad = np.pad(x, (pad_left, pad_right), mode="edge")
    kernel = np.ones(window) / window
    return np.convolve(xpad, kernel, mode="valid")


def compute_stride(universe, interval_ns):
    dt_ps = universe.trajectory.dt
    dt_ns = dt_ps / 1000.0
    stride = int(interval_ns / dt_ns) if dt_ns > 0 else 1
    return max(stride, 1), dt_ns


# Build the per-chain residue-pair list once, from the FIRST replica's
# universe (topology/residue numbering is identical across replicas).
u0 = mda.Universe(topology, trajfiles[0])
chain_res_pairs = {}   # chain_name -> list of (i, j) residue-index pairs
chain_n_res = {}
for chain_name in chain_groups:
    chain_atoms0 = u0.atoms[ndx[chain_name]]
    n_res = len(chain_atoms0.residues)
    pairs = [(i, j) for i, j in combinations(range(n_res), 2)
             if abs(i - j) > residue_exclude]
    chain_n_res[chain_name] = n_res
    chain_res_pairs[chain_name] = pairs
    print(f"\nChain '{chain_name}': {n_res} residues, {len(pairs)} intra-chain "
          f"residue pairs to evaluate per frame (RESIDUE_EXCLUDE={residue_exclude})")
    if len(pairs) == 0:
        print(f"  WARNING: no residue pairs survive RESIDUE_EXCLUDE={residue_exclude} "
              f"for chain '{chain_name}' (chain too short?) — its contribution will be 0.")

# =========================
# PLOT SETUP (single shared figure — one line per replica)
# =========================
plt.rcParams.update({"font.size": 16})
fig, ax = plt.subplots(figsize=(8, 6))

# =========================
# LOOP OVER REPLICAS
# =========================
for rep_name, tfile in zip(repnames, trajfiles):

    print(f"\nProcessing {rep_name}")
    u = mda.Universe(topology, tfile)

    chains = {}
    for chain_name in chain_groups:
        chain_atoms = u.atoms[ndx[chain_name]]
        residues = list(chain_atoms.residues)
        if len(residues) != chain_n_res[chain_name]:
            sys.exit(f"ERROR: residue count mismatch for chain '{chain_name}' between "
                      f"replicas ({rep_name} has {len(residues)}, expected "
                      f"{chain_n_res[chain_name]}) — topology should be identical "
                      f"across replicas.")
        chains[chain_name] = residues

    times = []
    total_contacts_per_frame = []

    outfile = os.path.join(outdir, f"{rep_name}_intra_contacts.dat")

    with open(outfile, "w") as f:
        header = "time_ns frame " + " ".join(chain_groups) + " total_contacts\n"
        f.write(header)

        stride, dt_ns = compute_stride(u, analysis_interval_ns)
        print(f"{rep_name} -> frame spacing: {dt_ns:.6f} ns | stride: {stride}")

        for ts in tqdm(u.trajectory[::stride], desc=rep_name):

            time_ns = ts.time / 1000.0
            times.append(time_ns)

            per_chain_contacts = {}
            total_contacts = 0

            for chain_name in chain_groups:
                residues = chains[chain_name]
                pairs = chain_res_pairs[chain_name]
                contacts = 0

                for i, j in pairs:
                    d = distances.distance_array(
                        residues[i].atoms.positions,
                        residues[j].atoms.positions,
                    )
                    if np.any(d < cutoff):
                        contacts += 1

                per_chain_contacts[chain_name] = contacts
                total_contacts += contacts

            total_contacts_per_frame.append(total_contacts)

            line = f"{time_ns:.4f} {ts.frame:d} "
            line += " ".join(f"{per_chain_contacts[c]:d}" for c in chain_groups)
            line += f" {total_contacts:d}\n"
            f.write(line)

    times = np.asarray(times)
    total_contacts_per_frame = np.asarray(total_contacts_per_frame)
    total_contacts_ma = moving_average(total_contacts_per_frame, ma_window)

    ax.plot(times, total_contacts_ma, label=rep_name, linewidth=2)

# =========================
# FINAL PLOT
# =========================
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Total Intra-chain Residue Contacts")
ax.legend()

plt.tight_layout()
plot_path = os.path.join(outdir, "intra_chain_contacts.png")
plt.savefig(plot_path, dpi=600)
plt.close(fig)

print("\nDone. Outputs saved in", outdir)