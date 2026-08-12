
CoCoMaps2 Figure 8 + Supplementary analysis pipeline
====================================================

Purpose
-------
Publication-oriented pipeline for the HPV16 E6/E7 peptide-HLA-A*02:01 CoCoMaps2 analysis.

Designed output:
Figure 8:
A. Interface area WT vs APL
B. Interaction network composition
C. Anchor position contribution (P1/P2/C-terminal)
D. HLA hotspot residue remodeling

Supplementary:
S1. Complete interaction type heatmap
S2. Peptide residue contact contribution
S3. HLA hotspot heatmap
S4. Pairwise WT-APL delta table
S5. Full interaction statistics table

Figure 7 color consistency
--------------------------

E6wt    #1A4F8A
E6apl   #8B3A0F

E7wt1   #5BA3D0
E7apl1  #E08A3C

E7wt2   #3A7EBF
E7apl2  #C1622A


Workflow
--------

1. Put CoCoMaps2 result zip files into:

raw_zip/

2. Run:

python3 scripts/01_extract_cocomaps2_results.py

3. Run:

python3 scripts/02_generate_figure8.py

4. Run:

python3 scripts/03_generate_supplementary_figures.py


Output
------

summary/
    cocomaps_master_summary.csv
    cocomaps_pairwise_delta.csv
    cocomaps_peptide_anchor.csv
    cocomaps_hla_hotspots.csv

figures/
    Figure8/
    Supplementary/


Recommended manuscript use
--------------------------

Main text:
Figure 8

Supplementary:
Detailed interaction tables and residue-level maps

Interpretation framework:
- E6apl: interaction remodeling
- E7apl1: enhanced anchoring but increased rigidity
- E7apl2: C-terminal anchor optimization
