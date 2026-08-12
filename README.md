# HPV16 E6/E7 Altered Peptide Ligand Design and Multi-Epitope CTL Validation

## Overview

This repository contains the analysis scripts, processed data, and figure-generation code associated with the study:

**“Immunoinformatics-Driven Rational Design of Altered Peptide Ligands and In Vitro Validation of a Multi-Epitope CTL Formulation Targeting HPV16 E6/E7”**

The study integrates immunoinformatics prediction, altered peptide ligand (APL) design, in vitro functional validation, molecular dynamics (MD) simulation, peptide–HLA interface analysis, HPV16 sequence-conservation analysis, and HLA population-coverage assessment.

The repository is organized to provide a direct correspondence between the datasets and scripts used to generate the main and supplementary figures reported in the manuscript.

---

## Repository Structure

```text
hpv16-immunoinformatics-apl-ctl/
│
├── README.md
├── LICENSE
├── requirements.txt
│
├── scripts/
│   ├── Figure1/
│   ├── Figure2/
│   ├── Figure3/
│   ├── Figure4/
│   ├── Figure5/
│   ├── Figure6/
│   │
│   ├── Figure7_MD/
│   │   ├── cpptraj/
│   │   ├── mmgbsa/
│   │   └── plotting/
│   │
│   ├── Figure8_COCOMAPS2/
│   │   ├── input_preparation/
│   │   ├── result_processing/
│   │   └── plotting/
│   │
│   └── Supplementary/
│       ├── FigureS1/
│       ├── FigureS2/
│       ├── FigureS3/
│       ├── FigureS4/
│       ├── FigureS5/
│
├── data/
│   ├── Figure1/
│   ├── Figure2/
│   ├── Figure3/
│   ├── Figure4/
│   ├── Figure5/
│   ├── Figure6/
│   │
│   ├── Figure7_MD/
│   │   ├── rmsd/
│   │   ├── rmsf/
│   │   ├── hbond/
│   │   └── mmgbsa/
│   │
│   ├── Figure8_COCOMAPS2/
│   │
│   └── Supplementary/
│       ├── HPV16_conservation/
│       ├── MD_block_analysis/
│       ├── IEDB_population_coverage/
│       └── Tables/
│
└── figures/
    ├── Figure1/
    ├── Figure2/
    ├── Figure3/
    ├── Figure4/
    ├── Figure5/
    ├── Figure6/
    ├── Figure7/
    ├── Figure8/
    └── Supplementary/
```

The `scripts/` directory contains analysis and plotting code, the `data/` directory contains the processed numerical data required to reproduce the figures, and the `figures/` directory contains representative final outputs.

---

# Study Design

The study focuses on HPV16 E6- and E7-derived HLA-A*02:01-restricted epitopes and rationally designed altered peptide ligands.

The wild-type peptide candidates analyzed in the study include:

| Peptide | Sequence   | HPV16 protein |
| ------- | ---------- | ------------- |
| E6wt    | TIHDIILECV | E6            |
| E7wt1   | LLMGTLGIV  | E7            |
| E7wt2   | YMLDLQPETT | E7            |

Selected altered peptide ligands include:

| Peptide | Sequence   |
| ------- | ---------- |
| E6apl   | YLHDIILECV |
| E7apl1  | LLMGTLGIY  |
| E7apl2  | YMLDLQPETV |
| E7apl3  | YLLDLQPETV |

The final multi-epitope formulation evaluated experimentally consisted of six peptides:

```text
E6wt
E7wt1
E7wt2
E6apl
E7apl2
E7apl3
```

E7apl1 was not included in the final formulation because improved predicted peptide–HLA binding was not accompanied by enhanced downstream effector activity.

---

# Main Figures

## Figure 1

Figure 1 summarizes the overall study design and computational-to-experimental workflow.

Associated materials are located in:

```text
scripts/Figure1/
data/Figure1/
figures/Figure1/
```

Where applicable, source data or figure-generation files used for the schematic are provided.

---

## Figure 2

Figure 2 presents the position-mapping analysis of wild-type HPV16 E6/E7 peptides and their corresponding altered peptide ligands (APLs). The figure illustrates the locations of the selected epitopes within the E6/E7 proteins and the residue substitutions introduced during APL design.

Associated materials are located in:

```text
scripts/Figure2/
data/Figure2/
figures/Figure2/
```

---

## Figure 3

Figure 3 contains experimental analyses related to peptide-induced CD8+ T-cell expansion and cell viability.

Processed quantitative data used for plotting are provided in:

```text
data/Figure3/
```

Figure-generation scripts are provided in:

```text
scripts/Figure3/
```

Representative microscopy images included in the publication may also be provided where appropriate. Large or personally identifiable raw experimental files are not included.

---

## Figure 4

Figure 4 contains IFN-γ ELISpot analyses evaluating peptide-specific T-cell responses.

The repository provides the numerical data used for statistical analysis and visualization, together with the corresponding plotting scripts.

```text
scripts/Figure4/
data/Figure4/
figures/Figure4/
```

Where applicable, background-subtracted ELISpot values and statistical-analysis inputs are included.

---

## Figure 5

Figure 5 contains peptide–HLA-A*02:01 stabilization and associated flow-cytometry analyses.

Processed MFI-based peptide-stabilization data and plotting scripts are available in:

```text
scripts/Figure5/
data/Figure5/
figures/Figure5/
```

Large raw flow-cytometry files are not necessarily included in the repository because the repository is intended primarily to reproduce the quantitative analyses and publication figures.

---

## Figure 6

Figure 6 contains cytotoxicity analyses of multi-epitope HPV-specific CTLs.

Processed donor-level cytotoxicity data and corresponding plotting scripts are provided in:

```text
scripts/Figure6/
data/Figure6/
figures/Figure6/
```

The published cytotoxicity analysis includes biological replication across independent donors where indicated in the manuscript.

---

# Figure 7 — Molecular Dynamics Analysis

Figure 7 summarizes the molecular dynamics analyses of peptide–HLA-A*02:01 complexes.

The six principal systems are:

```text
E6wt
E6apl
E7wt1
E7apl1
E7wt2
E7apl2
```

The MD analysis includes:

* root-mean-square deviation (RMSD);
* root-mean-square fluctuation (RMSF);
* peptide–HLA hydrogen-bond analysis;
* MM-GBSA binding-energy analysis.

Scripts are organized under:

```text
scripts/Figure7_MD/
├── cpptraj/
├── mmgbsa/
└── plotting/
```

Processed numerical data are provided under:

```text
data/Figure7_MD/
├── rmsd/
├── rmsf/
├── hbond/
└── mmgbsa/
```

## MD trajectory analysis

The equilibrated portion of each trajectory from **20 to 50 ns** was used for the analyses reported in the manuscript.

Trajectory processing and structural analyses were performed using `cpptraj`.

Relevant trajectory-processing scripts and analysis commands are included in:

```text
scripts/Figure7_MD/cpptraj/
```

## MM-GBSA analysis

Binding-energy calculations were performed using `MMPBSA.py`.

The generalized Born model and ionic-strength parameters used for the reported analyses were:

```text
igb = 5
saltcon = 0.15 M
```

Corresponding scripts and processed outputs are available in:

```text
scripts/Figure7_MD/mmgbsa/
data/Figure7_MD/mmgbsa/
```

## Raw trajectories

Full MD trajectory files are not included in this GitHub repository because of their large file sizes.

Instead, the processed numerical data required to reproduce the published plots are provided.

---

# Figure 8 — COCOMAPS 2.0 Interface Analysis

Figure 8 characterizes peptide–HLA-A*02:01 interface remodeling using **COCOMAPS 2.0**.

The analysis compares:

```text
E6wt vs E6apl
E7wt1 vs E7apl1
E7wt2 vs E7apl2
```

The analysis includes:

* peptide–HLA interface area;
* favorable interfacial interactions;
* anchor-residue interaction contributions;
* HLA groove hotspot contact patterns.

Scripts are organized under:

```text
scripts/Figure8_COCOMAPS2/
├── input_preparation/
├── result_processing/
└── plotting/
```

Processed COCOMAPS 2.0 summary data are stored in:

```text
data/Figure8_COCOMAPS2/
```

Key processed files include:

```text
cocomaps_master_summary.csv
cocomaps_anchor_contact.csv
cocomaps_hla_hotspots.csv
```

The final Figure 8 plotting script is located under:

```text
scripts/Figure8_COCOMAPS2/plotting/
```

COCOMAPS 2.0 input structures were generated from representative peptide–HLA-A*02:01 structures extracted from the MD analysis.

---

# Supplementary Analyses

## Figures S1 and S2 — HPV16 Sequence Conservation

HPV16 E6/E7 sequence conservation was evaluated across representative HPV16 lineages and sublineages.

The analysis included 16 representative HPV16 isolates spanning major viral lineages/sublineages.

The sequence-analysis workflow consisted of:

```text
sequence acquisition
        ↓
E6/E7 sequence extraction
        ↓
translation
        ↓
multiple-sequence alignment
        ↓
residue-level conservation analysis
        ↓
epitope and APL anchor-site analysis
```

Multiple-sequence alignment was performed using:

```text
MAFFT v7.453
```

The analysis pipeline is provided under:

```text
scripts/Supplementary/HPV16_conservation/
```

Processed conservation datasets are available in:

```text
data/Supplementary/HPV16_conservation/
```

Major output files include:

```text
full_length_conservation_E6.tsv
full_length_conservation_E7.tsv
epitope_conservation.tsv
APL_anchor_conservation.tsv
```

These files contain the numerical data used to generate Figures S1 and S2 and the associated supplementary conservation table.

---

## Figure S3 — COCOMAPS 2.0 Supplementary Analysis

Figure S3 provides supplementary visualization of peptide–HLA interface features derived from the COCOMAPS 2.0 analysis.

Associated plotting scripts are located in:

```text
scripts/Supplementary/FigureS3/
```

Because Figure S3 and Figure 8 originate from the same COCOMAPS 2.0 analysis pipeline, shared processed datasets are stored centrally under:

```text
data/Figure8_COCOMAPS2/
```

This avoids unnecessary duplication of identical source data.

---

## Figure S4 — MD Block-Averaging Analysis

To evaluate the temporal robustness of the molecular dynamics results, the equilibrated **20–50 ns** interval of each trajectory was divided into three consecutive 10 ns blocks:

```text
20–30 ns
30–40 ns
40–50 ns
```

RMSF and MM-GBSA analyses were independently performed for each block.

### RMSF block analysis

Residue-level backbone RMSF values were calculated independently for the three trajectory blocks.

The final plots report the block-averaged RMSF together with the standard error of the mean across blocks.

Processed RMSF summary files include:

```text
E6wt_10mer_rmsf_meanSEM.dat
E6apl_10mer_rmsf_meanSEM.dat

E7wt1_9mer_rmsf_meanSEM.dat
E7apl1_9mer_rmsf_meanSEM.dat

E7wt2_10mer_rmsf_meanSEM.dat
E7apl2_10mer_rmsf_meanSEM.dat
```

### MM-GBSA block analysis

MM-GBSA calculations were independently performed for each 10 ns trajectory block using the same parameters as the full-trajectory analysis.

Processed block-level MM-GBSA data and plotting scripts are located under:

```text
scripts/Supplementary/FigureS4/
data/Supplementary/MD_block_analysis/
```

The block analysis was designed as a temporal robustness assessment and should not be interpreted as a substitute for fully independent replicate MD simulations.

---

## Figure S5 — HLA Population Coverage

HLA population coverage was evaluated using the **IEDB Population Coverage Tool**.

The analysis was based on the peptide/HLA settings described in the manuscript.

Because the IEDB population-coverage analysis was performed using the web-based IEDB tool, the repository provides:

* analysis inputs;
* peptide/HLA settings;
* selected population datasets or regions;
* processed output tables;
* figure-generation data;
* documentation describing how the analysis was performed.

These materials are located under:

```text
data/Supplementary/IEDB_population_coverage/
scripts/Supplementary/FigureS5/
```

---

# Supplementary Tables

## Table S1

Table S1 summarizes HPV16 sequence-conservation results associated with the E6/E7 conservation analysis.

The underlying numerical data are available in:

```text
data/Supplementary/HPV16_conservation/
```

---

## Table S2

Table S2 provides an assay overview matrix documenting which computational and experimental evaluations were performed for individual peptide candidates.

The machine-readable source data are stored under:

```text
data/Supplementary/Tables/
```

---

## Table S3

Table S3 summarizes the IEDB HLA population-coverage analysis.

Underlying data and analysis settings are available in:

```text
data/Supplementary/IEDB_population_coverage/
```

# Software and Tools

The analyses reported in this repository used a combination of computational and experimental-analysis tools, including:

* NetMHCpan 4.1
* Amber / AmberTools
* `cpptraj`
* `MMPBSA.py`
* COCOMAPS 2.0
* MAFFT v7.453
* IEDB Population Coverage Tool
* Python
* R

Exact analysis parameters and software-specific commands are provided in the corresponding scripts and README files within individual analysis directories.

---

# Python Dependencies

Python plotting and data-processing scripts primarily rely on standard scientific Python packages such as:

```text
numpy
pandas
matplotlib
```

Additional dependencies, where required, are documented in the relevant script directories.

A `requirements.txt` file may be used to install the common Python environment:

```bash
pip install -r requirements.txt
```

---

# Reproducibility

For most figures, the repository contains two components required for reproduction:

1. the processed numerical data used in the published figure;
2. the corresponding plotting or analysis script.

For example:

```text
data/Figure8_COCOMAPS2/
        ↓
scripts/Figure8_COCOMAPS2/plotting/
        ↓
figures/Figure8/
```

Large simulation trajectories and other high-volume intermediate files are not included when they are unnecessary for reproducing the published numerical plots.

The repository therefore focuses on providing the processed analysis outputs and scripts directly associated with the results reported in the manuscript.

---

# Data Organization Principles

To reduce unnecessary duplication:

* shared COCOMAPS 2.0 data used by Figure 8 and Figure S3 are stored under `data/Figure8_COCOMAPS2/`;
* shared HPV16 conservation data used by Figures S1, S2, and Table S1 are stored under `data/Supplementary/HPV16_conservation/`;
* shared MD block-analysis data used by Figure S4 are stored under `data/Supplementary/MD_block_analysis/`;
* IEDB population-coverage data used by Figure S5 and Table S3 are stored under `data/Supplementary/IEDB_population_coverage/`.

---

# Notes on Raw Data

This repository is intended to provide the scripts and processed numerical datasets necessary to reproduce the analyses and figures reported in the manuscript.

Large raw files, including full molecular-dynamics trajectories, may be excluded because of repository file-size constraints.

Where raw experimental image or flow-cytometry files are not included, the processed numerical values used for the published analyses are provided where appropriate.

---

# Citation

If you use the scripts or processed datasets provided in this repository, please cite the associated publication:

> Dong D, Zhang X, Li B.
> **Immunoinformatics-Driven Rational Design of Altered Peptide Ligands and In Vitro Validation of a Multi-Epitope CTL Formulation Targeting HPV16 E6/E7.**
> *Vaccines*.
> Publication details will be updated after final publication.

Repository:

```text
https://github.com/dongdian1022/hpv16-immunoinformatics-apl-ctl
```

---

# Key Methodological References

The molecular-dynamics trajectory analyses used `cpptraj`:

> Roe DR, Cheatham TE III. PTRAJ and CPPTRAJ: Software for Processing and Analysis of Molecular Dynamics Trajectory Data. *J Chem Theory Comput.* 2013;9:3084–3095.

MM-GBSA calculations were performed using `MMPBSA.py`:

> Miller BR III, McGee TD Jr, Swails JM, Homeyer N, Gohlke H, Roitberg AE. MMPBSA.py: An Efficient Program for End-State Free Energy Calculations. *J Chem Theory Comput.* 2012;8:3314–3321.

Additional software and methodological references are provided in the associated manuscript.

---

# Contact

For questions regarding the repository or analysis workflow, please contact:

**Dian Dong**
College of Life Sciences, University of Chinese Academy of Sciences
BGI Research, BGI-Shenzhen
Email: [dongdian19@mails.ucas.ac.cn](mailto:dongdian19@mails.ucas.ac.cn)
dongdian@genomics.cn

# License

Please refer to the `LICENSE` file for terms of reuse of the code and repository materials.
