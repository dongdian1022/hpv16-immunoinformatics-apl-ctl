# Molecular Dynamics (MD) Scripts

This directory contains molecular dynamics simulation input files and analysis scripts.

## Software

- Amber26 / AmberTools26
- pmemd.cuda (GPU-accelerated MD engine)
- MMPBSA.py (MM-GBSA calculations)
- cpptraj (trajectory analysis)

## Force Field and Parameters

- Force field: ff19SB
- Water model: OPC
- GB model: igb = 5 (Onufriev GB-OBC II)
- Salt concentration: 0.15 M

## Contents

- System preparation scripts (tleap inputs)
- Energy minimization and equilibration protocols
- Production MD input files (50 ns trajectories)
- cpptraj scripts for:
  - RMSD analysis
  - RMSF per-residue fluctuation
  - Hydrogen bond occupancy
- MM-GBSA calculation scripts

## Output

- Backbone RMSD/RMSF plots (Figure 7A–F)
- Hydrogen bond occupancy analysis (Figure 7G–I)
- Binding free energy decomposition (Table 2)

## Notes

All simulations were performed using a single-trajectory MM-GBSA approach. Entropy contributions were not included.
