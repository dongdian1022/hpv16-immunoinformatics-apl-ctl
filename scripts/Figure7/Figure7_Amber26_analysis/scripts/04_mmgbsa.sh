#!/bin/bash
# =============================================================================
# 04_mmgbsa.sh
# MM-GBSA binding free energy on the 20-50 ns equilibrated (stripped) trajectory.
# Topologies are split with ante-MMPBSA.py; energies with MMPBSA.py (igb=5, 0.15 M).
# 结合自由能计算 (MM-GBSA, igb=5, 盐浓度 0.15 M)，取 20-50 ns 平衡段。
#
# NOTE on masks: REC_MASK follows the receptor/ligand split used to generate the
#   published Table 2 values (:1-376, as in the v2 production run).
#
# Frame convention: clean traj = 20-50 ns already; here we use its full range.
#   (If you feed the FULL 0-50 ns clean traj, set startframe=4001 instead.)
#
# Input : topology/*.prmtop , analysis/clean_traj/<name>_clean_20_50ns.nc
# Output: analysis/mmpbsa/data , analysis/mmpbsa/logs , analysis/mmpbsa/prmtops
# Usage : bash 04_mmgbsa.sh
# =============================================================================
set -euo pipefail

AMBER_SH="${AMBER_SH:-/path/to/amber/amber.sh}"
source "$AMBER_SH"

mkdir -p analysis/mmpbsa/data analysis/mmpbsa/logs analysis/mmpbsa/prmtops

REC_MASK=":1-376"

for prmtop in topology/*.prmtop; do
    name=$(basename "$prmtop" .prmtop)
    echo "=== $name ==="

    # 1. Split topology into complex / receptor / ligand (solvent stripped)
    if [ ! -f "analysis/mmpbsa/prmtops/${name}_com.prmtop" ]; then
        echo "  [1/2] splitting topology"
        ante-MMPBSA.py -p topology/${name}.prmtop \
            -c analysis/mmpbsa/prmtops/${name}_com.prmtop \
            -r analysis/mmpbsa/prmtops/${name}_rec.prmtop \
            -l analysis/mmpbsa/prmtops/${name}_lig.prmtop \
            -s ":WAT,Na+,Cl-" -n "$REC_MASK" > /dev/null 2>&1
    fi

    # 2. MM-GBSA control file
    cat <<EOF > analysis/mmpbsa/mmpbsa_${name}.in
MM-GBSA on 20-50 ns segment
&general
   keep_files=0, interval=10,
/
&gb
   igb=5, saltcon=0.15,
/
EOF

    echo "  [2/2] MMPBSA.py"
    MMPBSA.py -O -i analysis/mmpbsa/mmpbsa_${name}.in \
        -o  analysis/mmpbsa/data/${name}_FINAL_energy.dat \
        -cp analysis/mmpbsa/prmtops/${name}_com.prmtop \
        -rp analysis/mmpbsa/prmtops/${name}_rec.prmtop \
        -lp analysis/mmpbsa/prmtops/${name}_lig.prmtop \
        -y  analysis/clean_traj/${name}_clean_20_50ns.nc \
        > analysis/mmpbsa/logs/${name}_mmpbsa.log 2>&1

    if [ $? -eq 0 ]; then echo "  OK: $name"; else echo "  FAIL: $name (see logs)"; fi
    rm -f analysis/mmpbsa/mmpbsa_${name}.in
done
echo "MM-GBSA complete. See analysis/mmpbsa/data/"
