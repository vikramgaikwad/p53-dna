#!/usr/bin/env python3
"""
Compute the mean-square displacement (MSD) of water molecules
that stay in the hydration layer of a protein.
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.lib.distances import distance_array
import matplotlib.pyplot as plt
from tqdm import tqdm

# ============================================================
# --- USER SETTINGS ---
# ============================================================
topology   = "../../prod-100fs-1ns.tpr"      # topology file (.tpr or .gro)
trajectory = "../../prod-100fs-1ns.xtc"      # trajectory file (.xtc)
protein_sel = "protein and not name H*" # heavy atoms of protein
water_sel   = "resname TIP3 and name OW" # GROMACS water oxygen
dc1, dc2 = 0.30, 0.50                   # nm: inner & outer cutoff
threshold = 0.9                         # must stay ≥ 90% frames in shell
outfile = "hydration_msd.txt"           # MSD output
plotfile = "hydration_msd.png"          # MSD plot
# ============================================================

print("Loading trajectory ...")
u = mda.Universe(topology, trajectory)
protein = u.select_atoms(protein_sel)
waters = u.select_atoms(water_sel)
n_frames = len(u.trajectory)
print(f"Loaded {n_frames} frames, {len(waters)} water oxygens.")

# Align trajectory to protein backbone to remove drift
print("Aligning trajectory to protein backbone ...")
align.AlignTraj(u, u, select=protein_sel, in_memory=True).run()

# --- STEP 1: Determine which waters stay in the hydration shell ---
print("Identifying waters that stay in hydration shell ...")
stay_count = np.zeros(len(waters), dtype=int)

for ts in tqdm(u.trajectory, desc="Scanning trajectory"):
    # compute minimum distance between each water and protein
    dist = distance_array(waters.positions, protein.positions, box=u.dimensions)
    min_d = dist.min(axis=1)
    stay_count += (min_d >= dc1) & (min_d <= dc2)

frac_in_shell = stay_count / n_frames
staying_idx = np.where(frac_in_shell >= threshold)[0]
hydration_waters = waters[staying_idx]
print(f"Waters staying ≥{threshold*100:.0f}% of time in shell: {len(hydration_waters)}")

if len(hydration_waters) == 0:
    raise SystemExit("No waters satisfy the staying condition; adjust cutoffs or threshold.")

# --- STEP 2: Collect positions of these waters across all frames ---
print("Collecting positions of hydration waters ...")
positions = []
for ts in u.trajectory:
    positions.append(hydration_waters.positions.copy())  # already in nm for GROMACS
positions = np.array(positions)  # shape = (n_frames, n_waters, 3)
n_frames, n_w, _ = positions.shape

# --- STEP 3: Compute MSD(t) ---
print("Computing MSD ...")
msd = np.zeros(n_frames)
for dt in tqdm(range(1, n_frames), desc="Lag times"):
    disp = positions[dt:] - positions[:-dt]
    sqdisp = np.sum(disp**2, axis=2)    # (n_frames - dt, n_w)
    msd[dt] = np.mean(sqdisp)

# --- STEP 4: Time axis ---
timestep = u.trajectory.dt / 1000.0   # ps -> ns
time = np.arange(n_frames) * timestep

# --- STEP 5: Save and plot ---
np.savetxt(outfile, np.column_stack([time, msd]), header="time(ns)  MSD(nm^2)")
print(f"MSD data saved to {outfile}")

plt.figure(figsize=(6,4))
plt.plot(time, msd, lw=2)
plt.xlabel("Time lag (ns)")
plt.ylabel("MSD (nm$^2$)")
plt.title("Hydration-shell water MSD")
plt.tight_layout()
plt.savefig(plotfile, dpi=200)
plt.show()

print(f"MSD plot saved as {plotfile}")

