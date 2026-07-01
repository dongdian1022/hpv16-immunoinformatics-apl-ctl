# Processed Data

This directory contains processed and analysis-ready datasets used for statistical evaluation and figure generation.

## Contents

- Background-subtracted ELISpot SFU/10⁶ cell values
- Aggregated flow cytometry summary statistics
- Fold expansion values for CD8⁺ T-cell assays
- Normalized HLA stabilization indices (T2 assay)
- Summary tables for MM-GBSA outputs

## Processing Steps

Raw data were processed as follows:
- Background subtraction (where applicable)
- Replicate averaging (technical replicates collapsed per donor)
- Normalization to control conditions
- Calculation of derived metrics (e.g., fold expansion, MFI index)

## Usage

Processed datasets are directly used in:
- `/scripts_r` for statistical analysis
- Figure generation for the manuscript

These datasets correspond to all quantitative results shown in Figures 3–6.
