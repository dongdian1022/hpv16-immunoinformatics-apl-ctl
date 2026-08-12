#!/bin/bash
# =============================================================================
# 03_analyze_trajectory.sh
# Post-processing with cpptraj on the 20-50 ns equilibrated segment.
# For each system: strip solvent -> write clean trajectory -> RMSD, RMSF, H-bonds.
# RMSF is computed AFTER superposing on the receptor (HLA) CA atoms, so it
# reflects true local peptide flexibility (removes global tumbling).
# RMSF 在以受体 CA 主链叠合后计算，反映多肽真实局部柔性。
#
# Frame convention: dt=0.002 ps, ntwx=2500 -> 1 frame / 5 ps -> 50 ns = 10000 frames.
#   20-50 ns segment = frames 4001..last.
# Residues: receptor (HLA) = 1-375 ; peptide = 376.. (9-mer ->384, 10-mer ->385).
#
# Input : topology/*.prmtop , traj/<name>_md.nc
# Output: analysis/clean_traj , analysis/rmsd , analysis/rmsf , analysis/hbond
# Usage : bash 03_analyze_trajectory.sh
# =============================================================================
set -euo pipefail

AMBER_SH="${AMBER_SH:-/path/to/amber/amber.sh}"
source "$AMBER_SH"

mkdir -p analysis/clean_traj analysis/rmsd analysis/rmsf analysis/hbond

REC_MASK=":1-375"          # HLA receptor
START_FRAME=4001           # 20 ns (at 5 ps/frame)

for prmtop in topology/*.prmtop; do
    name=$(basename "$prmtop" .prmtop)
    echo ">> Analyzing $name (20-50 ns)"

    if [[ $name == *"10mer"* ]]; then
        PEP_MASK=":376-385"; MAX_RES="385"
    else
        PEP_MASK=":376-384"; MAX_RES="384"
    fi

    cat <<EOF > cpptraj_tmp.in
parm $prmtop
trajin traj/${name}_md.nc $START_FRAME last

# center + remove periodicity, then strip solvent/ions
autoimage
strip :WAT,Na+,Cl-
trajout analysis/clean_traj/${name}_clean_20_50ns.nc

# superpose on receptor CA -> removes global tumbling
rms fit ${REC_MASK}@CA

# backbone RMSD of the whole complex (reference = first retained frame)
rms Backbone_RMSD :1-${MAX_RES}@CA out analysis/rmsd/${name}_rmsd.dat time 0.005

# per-residue peptide flexibility
atomicfluct Peptide_RMSF $PEP_MASK out analysis/rmsf/${name}_rmsf.dat byres

# intra-complex hydrogen bonds
hbond Complex_Hbonds :1-${MAX_RES} \
      out    analysis/hbond/${name}_hbond_time.dat \
      avgout analysis/hbond/${name}_hbond_avg.dat
run
quit
EOF

    cpptraj -i cpptraj_tmp.in > analysis/clean_traj/${name}_cpptraj.log 2>&1
    echo "   done: $name"
done

rm -f cpptraj_tmp.in
echo "Analysis complete. See analysis/ subfolders."
