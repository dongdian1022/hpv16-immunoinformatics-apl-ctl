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

---

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

---

## Repository Structure
