#!/bin/bash
# =============================================================================
# 01_build_topology.sh
# Build AMBER topology/coordinate files for HLA-A*02:01 / peptide complexes.
# Force field: ff19SB (protein) + OPC (water); truncated-octahedron solvent box.
# 力场: ff19SB (蛋白) + OPC (水); 截角八面体水盒。
#
# Input : pdb_clean/*.pdb   (cleaned complex structures)
# Output: topology/<name>.prmtop , coords/<name>.inpcrd , logs/build_<name>.log
# Usage : bash 01_build_topology.sh
# =============================================================================
set -euo pipefail

# ---- User-configurable AMBER environment --------------------------------------
AMBER_SH="${AMBER_SH:-/path/to/amber/amber.sh}"   # edit or export AMBER_SH
source "$AMBER_SH"

mkdir -p logs topology coords scripts

shopt -s nullglob
pdb_list=(pdb_clean/*.pdb)
if [ ${#pdb_list[@]} -eq 0 ]; then
    echo "ERROR: no .pdb files found in pdb_clean/"; exit 1
fi

for pdb_file in "${pdb_list[@]}"; do
    name=$(basename "$pdb_file" .pdb)
    echo ">> Building system: $name"

    cat <<EOF > scripts/build_${name}.in
source leaprc.protein.ff19SB
source leaprc.water.opc
complex = loadpdb $pdb_file
addions complex Na+ 0
addions complex Cl- 0
solvateOct complex OPCBOX 10.0
saveamberparm complex topology/${name}.prmtop coords/${name}.inpcrd
quit
EOF

    tleap -s -f scripts/build_${name}.in > logs/build_${name}.log 2>&1

    if grep -q "Errors = 0" logs/build_${name}.log; then
        echo "   OK  ($name, Errors = 0)"
    else
        echo "   FAIL ($name) -> see logs/build_${name}.log"
    fi
done

echo "Done. Topologies in topology/ , coordinates in coords/"
