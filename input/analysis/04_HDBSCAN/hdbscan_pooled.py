#!/usr/bin/env python3
"""
Pooled HDBSCAN clustering across all replicas of the Aβ42 monomer pipeline.

Consumes every replica's production trajectory at once (not one job per
replica) and produces a single combined clustering result:
  - representative (medoid) PDB per cluster, written from the FULL topology
  - PCA scatter plot
  - labels.npy / pca_coords.npy

All configuration comes from environment variables (set by the calling
script_hdbscan_pooled.sh launcher) rather than a hardcoded USER SETTINGS
block, matching the env-var injection pattern used elsewhere in the pipeline.

Required env vars:
    TOPFILE     - path to nowater.gro (full-topology, water-stripped)
    NDXFILE     - path to index_analysis.ndx (must contain INDEX_GROUP_NAME)
    TRAJFILES   - comma-separated list of trajectory paths, one per replica

Optional env vars (defaults shown):
    OUTDIR                 = reps_monomer
    INDEX_GROUP_NAME       = Protein
    ATOM_NAME_FILTER       = CA            ("CA", "CA or CB", "all")
    PCA_COMPONENTS         = 10
    MIN_CLUSTER_FRACTION   = 0.001
    STRIDE                 = 1             (applied uniformly to all trajectories)
    FRAME_START            = 0             (applied uniformly, after stride)
    FRAME_END              = "" (empty = end of trajectory)
    SAVE_CLUSTER_FRAMES    = 0             (1 to export every frame of CLUSTERS_TO_EXPORT)
    CLUSTERS_TO_EXPORT     = 0             (comma-separated cluster ids)
"""

import sys
import os
import numpy as np
import mdtraj as md
from sklearn.decomposition import PCA
import hdbscan
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import Counter
from scipy.spatial.distance import cdist
from tqdm import tqdm


def env(name, default=None, cast=str):
    val = os.environ.get(name)
    if val is None or val == "":
        return default
    return cast(val)


def env_bool(name, default=False):
    val = os.environ.get(name)
    if val is None or val == "":
        return default
    return val.strip().lower() in ("1", "true", "yes", "on")


# ====================== CONFIG (from environment) ======================
topfile = env("TOPFILE")
ndxfile = env("NDXFILE")
trajfiles_raw = env("TRAJFILES")

if not topfile:
    sys.exit("ERROR: TOPFILE env var must be set (path to nowater.gro)")
if not ndxfile:
    sys.exit("ERROR: NDXFILE env var must be set (path to index_analysis.ndx)")
if not trajfiles_raw:
    sys.exit("ERROR: TRAJFILES env var must be set (comma-separated trajectory paths)")

trajfiles = [t.strip() for t in trajfiles_raw.split(",") if t.strip()]
if len(trajfiles) == 0:
    sys.exit("ERROR: TRAJFILES resolved to zero trajectories")

outdir = env("OUTDIR", "reps_monomer")
index_group_name = env("INDEX_GROUP_NAME", "Protein")
atom_name_filter = env("ATOM_NAME_FILTER", "CA")
pca_components = env("PCA_COMPONENTS", 10, int)
min_cluster_fraction = env("MIN_CLUSTER_FRACTION", 0.01, float)
stride = env("STRIDE", 1, int)
frame_start = env("FRAME_START", 0, int)
frame_end_raw = env("FRAME_END", "")
frame_end = int(frame_end_raw) if frame_end_raw else None
save_cluster_frames = env_bool("SAVE_CLUSTER_FRAMES", False)
clusters_to_export = [int(c) for c in env("CLUSTERS_TO_EXPORT", "0").split(",") if c.strip() != ""]

os.makedirs(outdir, exist_ok=True)

print("=== Pooled HDBSCAN — configuration ===")
print(f"  TOPFILE              : {topfile}")
print(f"  NDXFILE              : {ndxfile}")
print(f"  Trajectories ({len(trajfiles)})   :")
for t in trajfiles:
    print(f"    - {t}")
print(f"  OUTDIR               : {outdir}")
print(f"  INDEX_GROUP_NAME     : {index_group_name}")
print(f"  ATOM_NAME_FILTER     : {atom_name_filter}")
print(f"  PCA_COMPONENTS       : {pca_components}")
print(f"  MIN_CLUSTER_FRACTION : {min_cluster_fraction}")
print(f"  STRIDE               : {stride}")
print(f"  FRAME_START/END      : {frame_start} / {frame_end if frame_end is not None else '(end)'}")
print(f"  SAVE_CLUSTER_FRAMES  : {save_cluster_frames}  (clusters: {clusters_to_export})")

# ====================== ATOM SELECTION ======================
print("\n=== Parsing GROMACS Index & Selecting Atoms ===")


def parse_ndx(ndx_filename):
    """Custom native reader to parse GROMACS index (.ndx) files."""
    groups = {}
    current_group = None
    with open(ndx_filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('['):
                current_group = line.strip('[] ').strip()
                groups[current_group] = []
            elif current_group is not None:
                # Convert 1-based GROMACS indices to 0-based Python indices
                groups[current_group].extend([int(x) - 1 for x in line.split()])
    return groups


# Load topology once — reused for both atom-name filtering and PDB export
top = md.load_topology(topfile)

try:
    ndx_groups = parse_ndx(ndxfile)
except IOError:
    raise FileNotFoundError(f"Could not find or open the index file: {ndxfile}")

if index_group_name not in ndx_groups:
    raise KeyError(
        f"Group '{index_group_name}' not found in {ndxfile}. "
        f"Available groups: {list(ndx_groups.keys())}"
    )

base_idx = ndx_groups[index_group_name]

if atom_name_filter.lower() != "all":
    names_query = " or ".join([f"name {n.strip()}" for n in atom_name_filter.split("or")])
    filter_idx = top.select(names_query)
    atom_idx = np.intersect1d(base_idx, filter_idx)
    print(f"Intersected GROMACS group '{index_group_name}' with MDTraj filter '{names_query}'")
else:
    atom_idx = np.array(base_idx)

atom_idx = np.sort(atom_idx).tolist()
print(f"Final count: Selected {len(atom_idx)} atoms for clustering.")

if len(atom_idx) == 0:
    raise ValueError("ATOM SELECTION ERROR: 0 atoms selected. Check your index file and filter.")

# ====================== TRAJECTORY LOADING ======================
print("\n=== Loading Trajectories ===")
traj_id = []
traj_segments = []
frame_mapping = []   # (traj_file, disk_frame_index) per loaded frame

for tfile in trajfiles:
    print(f"\nProcessing : {tfile}")

    with md.open(tfile) as f:
        total_file_frames = len(f)

    start_idx = frame_start // stride
    end_disk = frame_end if frame_end is not None else total_file_frames
    end_idx = min(end_disk, total_file_frames) // stride
    n_expected = max(0, end_idx - start_idx)

    print(f"  Range      : {frame_start}:{end_disk if frame_end is not None else 'end'}  |  Stride: {stride}")

    with tqdm(total=n_expected, desc="  Loading frames", unit="fr", leave=True) as pbar:
        tr_full = md.load(tfile, top=topfile, atom_indices=atom_idx, stride=stride)
        tr2 = tr_full[start_idx:end_idx]
        pbar.update(tr2.n_frames)

    print(f"  → Loaded {tr2.n_frames} frames (stride={stride})")
    traj_segments.append(tr2)
    traj_id.extend([tfile] * tr2.n_frames)

    for i in range(start_idx, start_idx + tr2.n_frames):
        frame_mapping.append((tfile, i * stride))

subtraj = traj_segments[0].join(traj_segments[1:]) if len(traj_segments) > 1 else traj_segments[0]
print(f"\nTOTAL FRAMES LOADED: {subtraj.n_frames}")
traj_id = np.array(traj_id)

if subtraj.n_frames == 0:
    sys.exit("ERROR: 0 frames loaded across all trajectories — check FRAME_START/FRAME_END/STRIDE.")

# ====================== ALIGNMENT ======================
print("\nAligning trajectory to first frame...")
with tqdm(total=subtraj.n_frames, desc="  Superposing", unit="fr", leave=True) as pbar:
    subtraj.superpose(subtraj[0])
    pbar.update(subtraj.n_frames)

# ====================== FEATURES ======================
print("\n=== Building Feature Matrix ===")
X = subtraj.xyz.reshape(subtraj.n_frames, -1)
X = X - X.mean(axis=0)
print(f"  Feature matrix shape: {X.shape}")

# ====================== PCA ======================
print("\n=== Running PCA ===")
n_components = min(pca_components, X.shape[0], X.shape[1])
with tqdm(total=1, desc="  Fitting PCA", unit="step", leave=True) as pbar:
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X)
    pbar.update(1)
print(f"  Explained variance per PC: "
      f"{ {i+1: round(v*100, 2) for i, v in enumerate(pca.explained_variance_ratio_)} }")

# ====================== HDBSCAN ======================
print("\n=== Running HDBSCAN ===")
min_cluster_size = max(5, int(min_cluster_fraction * subtraj.n_frames))
print(f"  min_cluster_size = {min_cluster_size}")

with tqdm(total=1, desc="  Clustering", unit="step", leave=True) as pbar:
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        metric="euclidean",
        core_dist_n_jobs=-1
    )
    labels = clusterer.fit_predict(X_pca)
    pbar.update(1)

# Relabel clusters by decreasing size
counts = Counter(labels)
counts.pop(-1, None)
sorted_clusters = [lab for lab, _ in counts.most_common()]
relabel_map = {old: new for new, old in enumerate(sorted_clusters)}
labels = np.array([relabel_map.get(lab, -1) for lab in labels])

n_clusters = len(set(labels) - {-1})
n_noise = int(np.sum(labels == -1))
print(f"  Clusters found: {n_clusters}  |  Noise points: {n_noise} ({100*n_noise/len(labels):.1f}%)")
print(f"  Cluster populations: {Counter(labels)}")

# ====================== MEDOID REPRESENTATIVES ======================
print("\n=== Saving Representative Structures ===")


def medoid_index_pca(X_pca, idxs):
    """Find medoid in PCA space — faster than full coordinate space,
    and equivalent for ranking because PCA is a linear isometry up to scaling."""
    sub = X_pca[idxs]
    D = cdist(sub, sub, metric='euclidean')
    return idxs[np.argmin(D.sum(axis=0))]


cluster_ids = sorted(lab for lab in set(labels) if lab != -1)

for lab in tqdm(cluster_ids, desc="  Writing representative PDBs", unit="cluster", leave=True):
    idxs = np.where(labels == lab)[0]
    med = medoid_index_pca(X_pca, idxs)
    outpdb = os.path.join(outdir, f"cluster_{lab}_rep.pdb")

    target_tfile, target_frame = frame_mapping[med]
    full_frame = md.load_frame(target_tfile, index=target_frame, top=topfile)
    full_frame.save_pdb(outpdb)

    tqdm.write(
        f"  Cluster {lab:2d}: size={len(idxs):4d}  "
        f"medoid={med} (Frame {target_frame} of {os.path.basename(target_tfile)})"
        f"  → {outpdb}"
    )

# ====================== FULL-FRAME EXPORT (optional clusters) ======================
if save_cluster_frames:
    print(f"\n=== Exporting all frames for clusters {clusters_to_export} ===")

    for lab in clusters_to_export:
        if lab not in cluster_ids:
            print(f"  Cluster {lab} not found — skipping.")
            continue

        idxs = np.where(labels == lab)[0]
        frame_dir = os.path.join(outdir, f"cluster_{lab}_frames")
        os.makedirs(frame_dir, exist_ok=True)
        print(f"  Cluster {lab}: {len(idxs)} frames → {frame_dir}/")

        for rank, idx in enumerate(tqdm(idxs, desc=f"  Writing cluster {lab} frames",
                                        unit="fr", leave=True)):
            target_tfile, target_frame = frame_mapping[idx]
            outpdb = os.path.join(frame_dir, f"cluster_{lab}_frame{rank:05d}.pdb")
            full_frame = md.load_frame(target_tfile, index=target_frame, top=topfile)
            full_frame.save_pdb(outpdb)

        print(f"  → Done. {len(idxs)} PDBs written to {frame_dir}/")

# ====================== SAVE ARRAYS ======================
print("\nSaving PCA coordinates and labels...")
np.save("labels.npy", labels)
np.save("pca_coords.npy", X_pca)

# ====================== PCA SCATTER PLOT ======================
print("\n=== Generating PCA scatter plot ===")
plt.figure(figsize=(10, 8))
unique_labels = sorted(set(labels))
noise_mask = labels == -1
cluster_labels = [lab for lab in unique_labels if lab != -1]

cmap_name = "tab10" if len(cluster_labels) <= 10 else "tab20"
cmap = cm.get_cmap(cmap_name)
cluster_colors = {lab: cmap(i / max(1, len(cluster_labels) - 1))
                  for i, lab in enumerate(cluster_labels)}

if np.any(noise_mask):
    plt.scatter(
        X_pca[noise_mask, 0], X_pca[noise_mask, 1],
        s=20, alpha=0.4, color="#BBBBBB", label="Noise"
    )

for lab in sorted(cluster_labels, key=lambda l: np.sum(labels == l)):
    mask = labels == lab
    plt.scatter(
        X_pca[mask, 0], X_pca[mask, 1],
        s=40, facecolors=cluster_colors[lab],
        edgecolors="black", linewidth=0.5,
        alpha=0.9, label=f"Cluster {lab}"
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend()
plt.tight_layout()
plt.savefig("hdbscan_pca_scatter.png", dpi=300)
plt.close()

# ====================== FINAL SUMMARY ======================
print("\n=== DONE ===")
print(f"Representative PDBs saved to : {outdir}/")
if save_cluster_frames:
    for lab in clusters_to_export:
        print(f"All frames, cluster {lab}     : {outdir}/cluster_{lab}_frames/")
print(f"PCA scatter plot saved as     : hdbscan_pca_scatter.png")
print(f"Labels saved as               : labels.npy")
print(f"PCA coords saved as           : pca_coords.npy")

print("\nCluster distribution per trajectory:")
for lab in sorted(set(labels)):
    if lab == -1:
        continue
    mask = labels == lab
    traj_counts = Counter(traj_id[mask])
    print(f"\n  Cluster {lab}:")
    for t, c in traj_counts.items():
        print(f"    {t} : {c} frames")