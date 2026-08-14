# Aβ42 All-Atom MD Workflow (CHARMM36m, GROMACS)

End-to-end all-atom molecular dynamics workflow: starts from a PDB structure,
runs system preparation, energy minimisation, NVT/NPT equilibration and
production MD, then hands off to a separate analysis pipeline. This script works
for Aβ42, to change the system just use a different PDB structure from a 
different protein, however, you will need to change the index ChainA on run_setup.sh 
to read the correct number of atoms of the new system

Two scripts drive everything:

- **`run_setup.sh`** — master setup + simulation script (prep → min →
  NVT → NPT → production).
- **`run_setup_analysis.sh`** — master analysis script, run once production
  trajectories are ready.

The simulation stages (NPT + production) can run either as **SLURM jobs**
(the original, cluster-based workflow) or **locally**, sequentially, in your
own shell, with no SLURM involved. Both modes reuse the exact same
`grompp`/`mdrun` commands and the same segment-continuation logic — only the
submission/dependency mechanism differs.

---

## 1. Requirements

- GROMACS (tested with 2021.2; both CPU and GPU builds are supported)
- Bash
- For SLURM mode: a SLURM cluster with `sbatch`, and either a CPU or GPU
  partition set up (see [SLURM execution](#4-slurm-execution))
- The analysis scripts have their own dependencies — see the comments at the
  top of `run_setup_analysis.sh` and the scripts it calls.

---

## 2. Installation

Clone the repository and make sure the expected layout is present:

```
.
├── run_setup.sh
├── run_setup_analysis.sh
├── watch_npt_submit_prd.sh
├── mdp/
│   ├── ions_allatom.mdp
│   ├── min_allatom.mdp
│   ├── eq_nvt_allatom.mdp
│   ├── eq_npt_allatom.mdp
│   └── prd_allatom.mdp
└── input/
    ├── simulation/
    │   ├── min_allatom.sh
    │   ├── eq_NVT_allatom.sh
    │   ├── run_segments_local.sh
    │   ├── concat_trj_replica.sh
    │   ├── slurm_MD_eq_NPT_allatom.sh
    │   ├── slurm_MD_eq_NPT_allatom.dat
    │   ├── slurm_MD_prd_allatom.sh
    │   └── slurm_MD_prd_allatom.dat
    └── analysis/
        ├── 00_index/index_maker.sh
        ├── 01_Rg/ ...
        └── (one subfolder per analysis — see run_setup_analysis.sh)
```

Place your input PDB file in the repository root (default filename
`ab42_1z0q.pdb`, configurable — see below).

No build step is required; the scripts call your GROMACS installation
directly.

---

## 3. Configuration

All user-facing configuration lives in the **USER CONFIGURATION** block at
the top of `run_setup.sh`. Nothing below that block needs editing for
normal use.

### Core settings

| Variable | Meaning |
|---|---|
| `PDB` | Input PDB filename (must sit next to the script) |
| `SYS` | System name prefix used for output files |
| `GMX` | Path to the `gmx` executable used for local steps (prep, minimisation, NVT, and — in local mode — NPT/production) |
| `GMXLIB` | External force-field directory — see [Force fields and GMXLIB](#5-force-fields-and-gmxlib) |
| `BOX_PADDING`, `BOX_TYPE`, `ION_CONC` | Box/solvation/ionisation parameters |
| `NTHREADS` | Thread count for local steps |
| `MAXWARN` | `grompp -maxwarn` value |
| `FF`, `WM` | Force field / water model menu indices passed to `pdb2gmx` |
| `REPLICATES` | Number of independent production replicates (each gets its own `4_prd/repN/`) |

### Execution mode

```bash
EXECUTION_MODE=slurm   # or: local
```

- `slurm` (default) — original behaviour, unchanged: NPT and production run
  as self-requeuing SLURM jobs, with a background watcher chaining them.
- `local` — no SLURM: NPT then each production replicate run sequentially,
  in the foreground, in the current shell.

See [Local execution](#4-local-execution) and
[SLURM execution](#4-slurm-execution) below.

### Step counts and block counts

```bash
NPT_NSTEPS=50           # nsteps written into the NPT .mdp for this run
NPT_DT=0.002            # dt (ps) written into the NPT .mdp for this run
NPT_SEGMENTS=1          # number of NPT blocks (SLURM Segments= / local loop count)

PROD_NSTEPS=5000        # nsteps written into the production .mdp for this run
PROD_DT=0.002           # dt (ps) written into the production .mdp for this run
PROD_SEGMENTS=3         # number of production blocks per replicate
```

These control:

- The `nsteps` and `dt` values written into the **copies** of
  `eq_npt_allatom.mdp` / `prd_allatom.mdp` used for the run (the template
  files under `mdp/` are never modified).
- The number of NPT/production blocks — the SLURM `Segments=` value in
  `slurm_MD_eq_NPT_allatom.dat` / `slurm_MD_prd_allatom.dat`, or the number
  of sequential loop iterations in local mode.

One block/segment covers `nsteps × dt` of simulation time and continues
directly from the previous block's output structure — this is unchanged
in both execution modes. The defaults above reproduce the original
workflow's behaviour exactly; only change them if you want a different
segment length or a different number of blocks.

---

## 4. Local execution

Set `EXECUTION_MODE=local` and run:

```bash
bash run_setup.sh
```

Steps 0–5 (topology, box, solvation, ionisation, minimisation, NVT) are
identical to SLURM mode. For steps 6–7, the script calls
`input/simulation/run_segments_local.sh` directly — once for NPT, then once
per production replicate, one after another, entirely in the foreground.
Nothing is submitted to a queue and no background watcher is started.

`run_segments_local.sh` reimplements the same per-segment logic the SLURM
scripts use (segment numbering from existing `NNN.gro` files, `tinit`
computed as `nsteps × dt × (segment − 1)` read straight from the `.mdp`,
`grompp` + `mdrun` continuing from the previous segment) — just without
`sbatch`, scratch staging, or self-requeue. It is safe to re-run: any
segment whose `.gro` already exists is skipped.

Once production finishes, each replicate's trajectory is concatenated
automatically (`concat_trj_replica.sh`), producing
`4_prd/repN/traj_all_center_PBC.xtc` — the same file the analysis scripts
expect in SLURM mode.

Local mode is intended for workstations or interactive nodes; it uses the
`GMX`/`NTHREADS` values from the configuration block for every stage,
including NPT and production.

---

## 4. SLURM execution

Set `EXECUTION_MODE=slurm` (the default) and run:

```bash
./run_setup.sh
```

This is the original workflow, unchanged:

1. NPT equilibration is submitted via `slurm_MD_eq_NPT_allatom.sh`, which
   self-requeues until `Segments=` (from `slurm_MD_eq_NPT_allatom.dat`) is
   reached.
2. `watch_npt_submit_prd.sh` runs in the background, waits for NPT's last
   segment, then submits production for every replicate
   (`slurm_MD_prd_allatom.sh`, self-requeuing the same way).
3. Once each replicate's production finishes, the watcher concatenates its
   trajectory automatically.

The step-count/block-count configuration variables described above are
applied to the copies of the `.mdp`/`.dat` files placed in `3_npt/` and each
`4_prd/repN/` *before* the SLURM scripts are invoked, so the existing
dependency chain, self-requeue logic, and directory layout are completely
unaffected.

Monitor with:

```bash
tail -f watch_npt_submit_prd_*.log
squeue -u $USER
```

---

## 5. Force fields and GMXLIB

`GMXLIB` lets you point GROMACS at an external directory of force fields
without touching your installed GROMACS copy:

```bash
export GMXLIB="/path/to/my/gromacs-topologies"
```

This directory may contain additional or modified force-field folders
(`*.ff`), water models, residue definitions, or other GROMACS-compatible
topology files.

**Why use it:**

- Avoids modifying files inside the installed GROMACS distribution.
- Prevents custom force fields from being overwritten when GROMACS is
  updated or reinstalled.
- Keeps project-specific topology files separate from the system
  installation.
- Makes it easier to test, share, and document custom force-field setups.

**Notes:**

- `GMXLIB` should point to the directory *containing* your force-field
  folders, not necessarily directly to a single `.ff` folder — e.g. if you
  have `~/my-ffs/charmm36-custom.ff/`, set
  `GMXLIB=~/my-ffs`.
- Whichever force-field folder you select must have the internal structure
  GROMACS expects (`.ff` folder with `forcefield.itp`, etc.).
- You do not need to set `GMXLIB` at all if you're only using force fields
  already bundled with your GROMACS installation — just set `FF` to the
  corresponding menu index.
- This repository does **not** ship any external force fields. You are
  responsible for obtaining them from an appropriate source and complying
  with their license and citation requirements.
- `run_setup.sh` checks that `GMXLIB`, if set, actually exists, and
  will stop with a clear error if it doesn't. It also does a best-effort
  scan for available `.ff` directories and warns (without stopping) if
  `FF` doesn't appear to match anything found — the definitive check is
  always `pdb2gmx` itself at STEP 0, which will report its own clear error
  if the selected force field is invalid.

---

## 6. Analysis

Once production trajectories exist (`4_prd/repN/traj_all_center_PBC.xtc`),
run:

```bash
./run_setup_analysis.sh
```

See the comments at the top of that script for the individual analyses it
runs and the scripts it calls.

`run_setup_analysis.sh` has its own `EXECUTION_MODE`, independent of
`run_setup.sh`:

- `EXECUTION_MODE=ssh` (default) — dispatches analyses across a pool of
  machines (`MACHINES=(...)`) over SSH, `MAX_JOBS_PER_MACHINE` at a time per
  machine. Requires passwordless SSH to each machine.
- `EXECUTION_MODE=local` — no SSH, no machine pool. Everything runs on the
  current machine instead, `LOCAL_MAX_JOBS` analyses at a time. Set
  `LOCAL_MAX_JOBS=1` to run them strictly one after another, or raise it to
  run a few at once — useful if you only have your own workstation.

Each analysis produces a .png plot alongside its data — per-replica plots land in 4_prd/repN/analysis/<subfolder>/, and pooled analyses (e.g. HDBSCAN, interactions) land in 4_prd/Analysis/<subfolder>/

> **Note on `MACHINES`:** the default hostnames (`bio209`, `bio210`, ...)
> are specific to this repository's original cluster and won't exist
> elsewhere. Replace them with your own reachable hostnames, or just use
> `EXECUTION_MODE=local`.

---
