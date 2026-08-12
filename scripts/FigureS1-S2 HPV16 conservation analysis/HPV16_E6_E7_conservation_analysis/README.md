HPV16 E6/E7 conservation analysis v2.5 QC-fixed

Fixes from v2.4:
1. GenBank annotation fallback:
   - gene
   - product
   - note
   - protein_id

2. Missing lineage detection and QC report.

3. Pipeline error handling:
   - set -e
   - input/output validation

4. Supports 16 HPV16 sublineages:
   A1-A4, B1-B4, C1-C4, D1-D4

Workflow:
download -> translation -> MAFFT -> conservation analysis