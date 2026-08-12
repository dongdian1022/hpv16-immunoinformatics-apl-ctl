# Molecular dynamics & MM-GBSA scripts (Figure 7)

Scripts to reproduce the atomistic molecular dynamics (MD) simulations and
MM-GBSA binding free energy analysis of HLA-A\*02:01 / peptide complexes for the
three simulated wild-type / altered-peptide-ligand (WT/APL) pairs
(**E6wt/E6apl, E7wt1/E7apl1, E7wt2/E7apl2**) reported in Figure 7 and Table 2.

## Pipeline

| Step | Script | Tool | Purpose |
|------|--------|------|---------|
| 1 | `scripts/01_build_topology.sh` | tleap | Build solvated topology/coordinates (ff19SB + OPC, truncated octahedron) |
| 2 | `scripts/02_run_md_50ns.sh` | pmemd.cuda | Minimize → heat → equilibrate → 50 ns production MD |
| 3 | `scripts/03_analyze_trajectory.sh` | cpptraj | Strip solvent, RMSD, receptor-aligned peptide RMSF, H-bonds (20–50 ns) |
| 4 | `scripts/04_mmgbsa.sh` | ante-MMPBSA.py / MMPBSA.py | MM-GBSA (igb=5, 0.15 M) on the 20–50 ns segment |
| 5 | `scripts/05_extract_mmgbsa.sh` | awk | Summarize energies into a table / CSV |

Run in order from a working directory that contains `pdb_clean/` (cleaned input
structures). Each step creates the folders the next step needs.

## Software

- AMBER / AmberTools (ff19SB, OPC water, `pmemd.cuda`, `cpptraj`, MMPBSA.py)
- Bash; an NVIDIA GPU with CUDA for the MD stage
- Edit the `AMBER_SH` variable (or `export AMBER_SH=/path/to/amber.sh`) at the
  top of each script to point at your AMBER installation.

## Simulation protocol (summary)

- Force field ff19SB; OPC water in a 10 Å truncated-octahedron box; system
  neutralized with Na⁺/Cl⁻.
- Minimization (5000 steps) → heating 0→300 K over 100 ps (Langevin, γ=2 ps⁻¹)
  → 1 ns NPT equilibration → 50 ns NPT production (2 fs step, SHAKE, 8 Å cutoff,
  Monte-Carlo barostat), coordinates saved every 5 ps (10 000 frames).
- Analysis and MM-GBSA use the **20–50 ns equilibrated segment** (frames
  4001–10000). RMSF is computed after superposition on the receptor (HLA) CA
  atoms.

## Residue conventions

Receptor (HLA) = residues 1–375; peptide = residue 376 onward (9-mer → 384,
10-mer → 385). The MM-GBSA receptor/ligand split uses `:1-376` to match the
production run that generated the published Table 2 values.

## Data

Input structures, production trajectories, and raw analysis outputs are archived
on Zenodo (DOI: _to be inserted after publishing the deposit_).

## Citation

If you use these scripts, please cite the associated article (Vaccines, MDPI)
and the Zenodo archive.

## License

MIT (see `LICENSE`).
