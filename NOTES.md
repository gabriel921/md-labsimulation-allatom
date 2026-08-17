# A Student's Guide to This Aβ42 Molecular Dynamics Workflow

This document walks through *why* each step in this repository exists —
the chemistry and physics behind it — and *what* each analysis actually
tells you about the protein. It assumes you know basic molecular biology
and a little physical chemistry, but nothing about MD simulation.

---

## 1. What are we actually simulating, and why?

**Amyloid-β 42 (Aβ42)** is a 42-residue peptide fragment cleaved from a
larger membrane protein (APP, amyloid precursor protein). It is one of the
main protein species found in the amyloid plaques characteristic of
Alzheimer's disease.

The key chemical fact that makes Aβ42 interesting — and hard to study — is
that it is an **intrinsically disordered protein (IDP)**. Unlike a folded
enzyme with one well-defined 3D structure, a monomer of Aβ42 in solution
doesn't sit in a single conformation. It constantly interconverts between
many different shapes: extended coil, transient helices, transient
β-hairpins. This structural promiscuity is thought to be exactly what lets
it aggregate — some of those transiently sampled conformations are
"aggregation-prone," exposing hydrophobic patches that stick to the same
patches on another Aβ42 molecule, seeding the growth of oligomers and
eventually the β-sheet-rich fibrils found in plaques.

Two features of the sequence matter a lot for its behaviour:

- **The central hydrophobic core (CHC), residues 17–21 (LVFFA)** — a run
  of hydrophobic and aromatic side chains that plays a central role in
  self-association and in nucleating β-sheet structure.
- **The hydrophobic C-terminus (roughly residues 29–42)** — Aβ42 differs
  from the shorter, less aggregation-prone Aβ40 by two extra hydrophobic
  residues (Ile41, Ala42) at this end, which is a major reason Aβ42
  aggregates faster and more readily than Aβ40.

Because no single crystal structure represents an IDP, we can't just look
at "the structure" the way you would for a folded enzyme. Instead, we let
physics run: place the peptide in explicit water, integrate Newton's
equations of motion for a long time, and see what conformations it visits
and how often. That's what this whole repository does.

---

## 2. The physics: what an MD engine (GROMACS) is doing

At its core, molecular dynamics is very simple: for every atom, compute
the force on it (from bonded terms — bond lengths, angles, dihedrals —
and non-bonded terms — van der Waals and electrostatics), then use
Newton's second law (F = ma) to figure out how its position and velocity
change over a tiny time step (`dt`, typically 2 fs = 0.002 ps here).
Repeat millions of times.

The **force field** (CHARMM36m, used throughout this repository) is the
set of equations and parameters that define those bonded/non-bonded terms
for every atom type — essentially an empirical model of the potential
energy surface, fit to reproduce experimental and quantum-chemical data.
CHARMM36m specifically was reparametrized to fix known issues older force
fields had with disordered proteins (they tended to make IDPs too
compact/helical) — which is exactly why it's the right choice for Aβ42.

We use **explicit solvent** (real water molecules, not an implicit
continuum) because water's specific hydrogen-bonding behaviour around
hydrophobic and charged residues is central to how an IDP behaves — an
implicit solvent model would miss the very hydration effects that drive
the CHC's hydrophobic collapse.

---

## 3. Walking through the simulation stages

### Step 0 — Topology generation (`pdb2gmx`)

Takes the raw PDB structure and (a) assigns force-field parameters to
every atom, (b) adds hydrogens (the PDB file typically doesn't resolve
them), and (c) chooses protonation states for titratable residues (His,
Glu, Asp, etc.) at the force field's assumed pH. This produces the
**topology** (`topol.top`) — the full list of atoms, bonds, angles,
dihedrals, and charges that everything downstream depends on.

### Step 1 — Box building (`editconf`)

Defines a simulation box around the peptide with enough padding
(`BOX_PADDING`) that the peptide never "sees" its own periodic image
through the box wall — periodic boundary conditions mean the box repeats
infinitely in every direction, so if the padding is too small the peptide
would spuriously interact with copies of itself.

### Step 2 — Solvation (`solvate`)

Fills the box with explicit water molecules (SPC/E model here) around the
peptide.

### Step 3 — Ionisation (`genion`)

Adds Na⁺/Cl⁻ ions to (a) neutralise the system's net charge — required
for the electrostatics algorithm (PME) to work correctly under periodic
boundary conditions — and (b) bring the ionic strength up to a
physiologically realistic concentration (`ION_CONC`, default 0.1 M),
since screening of Aβ42's charged side chains by ions affects its
conformational ensemble.

### Step 4 — Energy minimisation

Before doing any dynamics, we need to relax the system — solvation and
ion placement create some atoms in awkward, high-energy positions (too
close to a neighbour, a slightly strained bond). Minimisation follows the
energy gradient downhill (steepest descent) to remove these clashes.
Skipping this step risks the simulation crashing in the first
femtoseconds from enormous forces.

### Step 5 — NVT equilibration

Now we start real dynamics, but under constant **N**umber of particles,
**V**olume, and **T**emperature — a thermostat is switched on and the
system is allowed to heat up to the target temperature and let velocities
equilibrate to a proper Maxwell-Boltzmann distribution, while the box
volume stays fixed. This step is short and mainly about "waking up" the
system's kinetic energy sensibly.

### Step 6 — NPT equilibration

Now volume is also allowed to fluctuate under a **barostat**, targeting a
constant pressure — this lets the box find its correct equilibrium
density (water's density at the given temperature/pressure), which the
fixed box from Step 1 doesn't automatically have. Only once density has
equilibrated is the system truly ready for production.

### Step 7 — Production MD

The actual "data-generating" run — long trajectories under NPT
(or NVT, depending on your `.mdp` settings), from which every analysis
in Section 6 below is computed.

**Segments/blocks:** production (and NPT) run in numbered chunks
(`001`, `002`, ...) rather than one single continuous run. Each block
continues exactly where the previous one left off (same coordinates,
`tinit` shifted forward by the elapsed time) — this exists purely so a
long simulation can be split across separate SLURM job submissions (each
with its own walltime limit) without losing continuity. Scientifically,
segmented and unsegmented runs of the same total length are equivalent.

**Replicates:** several independent copies of production
(`REPLICATES`) are run from the *same* NPT endpoint but will diverge
immediately, because SLURM/GROMACS reseeds random velocities each time —
this gives independent samples of the conformational ensemble, which
matters a lot for an IDP: a single trajectory, however long, can get
"stuck" for a while in one region of conformational space, so multiple
independent replicates are the standard way to check that your sampling
isn't an artefact of one particular trajectory's random walk.

---

## 4. The analyses — what each one tells you about the protein

Once you have trajectories, the raw coordinates aren't very interpretable
on their own — thousands of frames of xyz positions. Each analysis below
reduces that firehose of coordinates into a chemically meaningful number
or plot.

### Rg — Radius of gyration

A single number per frame: the mass-weighted RMS distance of every atom
from the molecule's centre of mass. Small Rg = compact/globular
conformation; large Rg = extended conformation. Tracking Rg over the
trajectory (and its distribution across replicates) tells you how
compact vs. extended Aβ42 tends to be, and whether it interconverts
between compact and extended states or stays in one regime.

### SASA — Solvent-accessible surface area

The surface area of the peptide that's in contact with (accessible to)
solvent. This directly reflects how much of the hydrophobic surface is
buried vs. exposed — a drop in SASA (especially of hydrophobic residues
specifically) usually signals hydrophobic collapse, e.g. the CHC burying
itself against the C-terminal hydrophobic stretch.

### DSSP — Secondary structure assignment

Classifies each residue, frame by frame, into a secondary-structure
category (α-helix, β-sheet/β-bridge, turn, bend, coil, ...) based on the
pattern of backbone hydrogen bonds. For an IDP like Aβ42, you're not
expecting a stable fold — you're looking at the *propensity* for
transient secondary structure at each residue, averaged over the
ensemble: which stretches spend more time as β-sheet/turn versus
featureless coil. This is one of the most direct read-outs of
aggregation-relevant conformational bias, since transient β-structure
(particularly around the CHC) is thought to be an early step toward
amyloid formation.

### HDBSCAN — Conformational clustering (pooled)

An unsupervised clustering algorithm applied across *all* replicates'
frames together (hence "pooled" — it needs the full combined ensemble to
find meaningful clusters). It groups structurally similar conformations
together and picks out the most heavily populated ones. Instead of
staring at thousands of individual frames, you get a handful of
representative structures with their relative populations — a practical
summary of "what does the conformational ensemble actually look like."

### Interactions / intra-contacts — Residue-residue contact maps (pooled)

Computes, across the ensemble, how often each pair of residues comes
into contact (within some distance cutoff) intramolecularly. Plotted as
a 2D contact map, this reveals which parts of the sequence tend to fold
back onto which other parts — e.g. persistent long-range contacts
between the CHC and the C-terminus would show up as an off-diagonal hot
spot, directly supporting (or not) a specific proposed collapse
mechanism.

---

## 5. Putting it together

None of these analyses stand alone — the real picture comes from
combining them. A representative "aggregation-prone-looking" snapshot
might show: low Rg (compact), reduced SASA at the CHC and C-terminus
specifically (buried against each other), a DSSP trace showing
above-baseline β/turn propensity in and around the CHC, a contact map
confirming CHC–C-terminus contacts, and that same conformational family
showing up as a well-populated HDBSCAN cluster. Building that kind of
composite argument from independent lines of evidence — rather than
trusting any single metric — is standard practice in IDP simulation work,
and is exactly the workflow this repository is set up to produce.

For a worked example of this kind of analysis on Aβ42 (and how
temperature modulates these same hydrophobic-collapse and β-sheet-forming
behaviours), see the published study linked in the main `README.md`.
