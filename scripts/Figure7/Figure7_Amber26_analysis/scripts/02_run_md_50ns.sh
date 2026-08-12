#!/bin/bash
# =============================================================================
# 02_run_md_50ns.sh
# Explicit-solvent MD for each HLA/peptide complex with pmemd.cuda (AMBER).
# Stages: minimization -> heating (0->300 K, 100 ps) -> NPT equilibration (1 ns)
#         -> production NPT (50 ns, 2 fs step, coords every 5 ps).
# 阶段: 能量最小化 -> 升温 -> 密度平衡 -> 50 ns 成品模拟。
#
# Input : topology/<name>.prmtop , coords/<name>.inpcrd
# Output: md_out/ , rst/ , traj/<name>_md.nc
# Usage : bash 02_run_md_50ns.sh              # runs all systems in topology/
#         SYSTEMS="E7wt1 E7apl1" bash 02_run_md_50ns.sh   # run a subset
# =============================================================================
set -euo pipefail

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"
mkdir -p scripts md_out rst traj

# ---- Shared MD control files --------------------------------------------------
cat <<'EOF' > scripts/min.in
Energy minimization
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=1, cut=8.0, ntpr=100,
 /
EOF

cat <<'EOF' > scripts/heat.in
Heating 0 -> 300 K (100 ps)
 &cntrl
  imin=0, irest=0, ntx=1,
  nstlim=50000, dt=0.002,
  ntc=2, ntf=2, cut=8.0, ntb=1,
  ntpr=1000, ntwx=1000,
  ntt=3, gamma_ln=2.0, tempi=0.0, temp0=300.0, ig=-1,
 /
EOF

cat <<'EOF' > scripts/eq.in
Density equilibration NPT (1 ns)
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=500000, dt=0.002,
  ntc=2, ntf=2, cut=8.0, ntb=2, ntp=1, barostat=2,
  ntpr=1000, ntwx=1000,
  ntt=3, gamma_ln=2.0, temp0=300.0, ig=-1,
 /
EOF

cat <<'EOF' > scripts/md.in
Production MD NPT (50 ns), coords every 5 ps
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=25000000, dt=0.002,
  ntc=2, ntf=2, cut=8.0, ntb=2, ntp=1, barostat=2,
  ntpr=2500, ntwx=2500,
  ntt=3, gamma_ln=2.0, temp0=300.0, ig=-1,
 /
EOF

# ---- System list --------------------------------------------------------------
if [ -n "${SYSTEMS:-}" ]; then
    read -r -a names <<< "$SYSTEMS"
else
    shopt -s nullglob
    names=(); for p in topology/*.prmtop; do names+=("$(basename "$p" .prmtop)"); done
fi

for name in "${names[@]}"; do
    if [ ! -f "topology/${name}.prmtop" ]; then
        echo "WARN: topology/${name}.prmtop not found, skipping."; continue
    fi
    echo "=== $name === $(date)"

    echo "[1/4] minimization"
    pmemd.cuda -O -i scripts/min.in  -p topology/${name}.prmtop -c coords/${name}.inpcrd \
        -o md_out/${name}_min.out  -r rst/${name}_min.rst  -ref coords/${name}.inpcrd

    echo "[2/4] heating"
    pmemd.cuda -O -i scripts/heat.in -p topology/${name}.prmtop -c rst/${name}_min.rst \
        -o md_out/${name}_heat.out -r rst/${name}_heat.rst -x traj/${name}_heat.nc -ref rst/${name}_min.rst

    echo "[3/4] equilibration"
    pmemd.cuda -O -i scripts/eq.in   -p topology/${name}.prmtop -c rst/${name}_heat.rst \
        -o md_out/${name}_eq.out   -r rst/${name}_eq.rst   -x traj/${name}_eq.nc

    echo "[4/4] production (50 ns)"
    pmemd.cuda -O -i scripts/md.in   -p topology/${name}.prmtop -c rst/${name}_eq.rst \
        -o md_out/${name}_md.out   -r rst/${name}_md.rst   -x traj/${name}_md.nc

    echo "OK: $name finished $(date)"
done
echo "All production runs complete."
