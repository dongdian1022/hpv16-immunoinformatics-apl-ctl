# hpv16-immunoinformatics-apl-ctl
Immunoinformatics-driven design and validation of HPV16 E6/E7 altered peptide ligands and multi-epitope CTL responses
# Immunoinformatics-Guided Design of Altered Peptide Ligands Targeting HPV16 E6/E7 and Experimental Validation of Multi-Epitope CTL Responses

## Overview
This repository contains all analysis scripts, processed datasets, and computational workflows associated with the study:
> *Immunoinformatics-driven rational design of altered peptide ligands (APLs) targeting HPV16 E6/E7 oncoproteins and experimental validation of multi-epitope cytotoxic T lymphocyte (CTL) responses.*
The study integrates:
- Immunoinformatics-based epitope prediction and prioritization
- Anchor-directed altered peptide ligand (APL) design
- Molecular dynamics (MD) simulations and MM-GBSA binding analysis
- Experimental validation using T2 binding assays, ELISpot, and CTL cytotoxicity assays

## Study Design Overview
The workflow is organized into two major phases:

### Phase I: In silico Epitope Identification and APL Design
- Sliding-window epitope scanning of HPV16 E6 and E7 proteins
- HLA-A*02:01 binding prediction using NetMHCpan-4.1
- Functional filtering using:
  - IEDB epitope database evidence
  - Literature-supported immunogenicity
  - Toxicity prediction (ToxinPred)
  - Antigenicity scoring (VaxiJen 3.0)
- Multi-anchor optimization at:
  - P1 (A-pocket)
  - P2 (B-pocket)
  - C-terminal residue (F-pocket)

### Phase II: Experimental and Structural Validation
- Peptide–HLA-A*02:01 binding (T2 stabilization assay)
- CD8⁺ T cell expansion assays
- IFN-γ ELISpot functional assays
- CTL cytotoxicity assays against HPV16 E7-expressing tumor cells
- Molecular dynamics simulations (Amber26) and MM-GBSA analysis

## Repository Structure

##################################################################################################################
HPV16 E6/E7 Immunoinformatics-Driven APL Design and Functional Validation Pipeline

This repository contains all data, analysis scripts, and structural modeling workflows associated with the study:

“Immunoinformatics-Based Rational Design of High-Affinity Altered Peptide Ligands and In Vitro Validation of a Multi-Epitope CTL Formulation Targeting HPV16 E6/E7”

The project integrates immunoinformatics epitope screening, anchor-optimized peptide design, in vitro functional immunology assays, and molecular dynamics (MD) simulations to construct and validate a multi-epitope CTL formulation targeting HPV16 E6/E7 oncoproteins.

1. Repository Structure
.
├── data_raw/              # Raw experimental and computational outputs
├── data_processed/        # Processed datasets used for figures and statistics
├── scripts_r/            # R scripts for statistical analysis and plotting (Fig 3–6)
├── scripts_md/           # Molecular dynamics setup, analysis, and MM-GBSA scripts (Fig 7)
├── figures_source/       # Source files for all figures (graph objects, panels, layouts)
├── results_tables/       # Final formatted tables used in manuscript (Table 1–2)
└── README.md
2. Data Overview
2.1 Experimental Data
T-cell expansion assays (3 stimulation rounds)
IFN-γ ELISpot assay (spot-forming units per 10⁶ cells)
T2 HLA-A*02:01 stabilization assay (MFI-based binding index)
Cytotoxicity assays (CFSE/7-AAD, multiple E:T ratios)
2.2 Computational Data
NetMHCpan-4.1 binding predictions
VaxiJen antigenicity scores
ToxinPred toxicity screening
Molecular dynamics trajectories (50 ns per system)
MM-GBSA binding free energy decomposition
3. Reproducibility Summary

All analyses were performed using:

R v4.4.3 (ggplot2 v4.0.3, dplyr v1.2.1)
AmberTools26 / AMBER 2026 (pmemd.cuda)
NetMHCpan-4.1
Python-based preprocessing utilities (custom scripts included)

Statistical analyses correspond exactly to those reported in Section 2.7 of the manuscript.

4. Key Outputs
4.1 Figures
Figure 3: CD8⁺ T-cell expansion dynamics
Figure 4: IFN-γ ELISpot functional responses
Figure 5: HLA-A*02:01 stabilization (T2 assay)
Figure 6: CTL cytotoxicity and specificity
Figure 7: MD simulations and binding mechanism analysis
4.2 Tables
Table 1: Epitope selection and immunoinformatics scoring
Table 2: MM-GBSA binding free energy decomposition
5. Figure-to-Script Mapping

Each figure is fully reproducible using the corresponding scripts:

Figure	Script	Description
Fig. 3	scripts_r/fig3_expansion.R	T-cell expansion analysis
Fig. 4	scripts_r/fig4_elispot.R	ELISpot SFU quantification
Fig. 5	scripts_r/fig5_t2_binding.R	HLA-A2 stabilization analysis
Fig. 6	scripts_r/fig6_cytotoxicity.R	Cytotoxicity and statistical testing
Fig. 7	scripts_md/	MD trajectories + RMSD/RMSF + MM-GBSA
6. Molecular Dynamics Pipeline

MD simulations were performed using:

HLA-A*02:01 templates (PDB: 1DUZ, 1I4F)
ff19SB force field
OPC water model
pmemd.cuda engine
50 ns production runs
20–50 ns equilibrium window used for analysis

Analysis includes:

Backbone RMSD
Per-residue RMSF
Hydrogen bond occupancy
MM-GBSA binding free energy (igb = 5, GB-OBC II)
7. Data Availability

All raw data, processed datasets, and analysis scripts will be:

deposited in a public repository (GitHub + Zenodo) upon publication
assigned a DOI via Zenodo for long-term archiving
8. Software Requirements
R environment
R ≥ 4.4.3
ggplot2 ≥ 4.0.3
dplyr ≥ 1.2.1
MD environment
AMBER 2026 / AmberTools26
CUDA-enabled pmemd.cuda recommended
9. Citation

If you use this repository, please cite the associated publication (DOI will be added upon publication).

10. License

To be specified upon publication (default: CC BY 4.0 recommended for Zenodo deposition).

11. Contact

For questions regarding data or code:

Corresponding author: dongdian@genomics.cn 
