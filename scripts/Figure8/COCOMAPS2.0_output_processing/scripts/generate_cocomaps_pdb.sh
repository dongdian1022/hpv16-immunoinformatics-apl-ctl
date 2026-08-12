#!/bin/bash

source /ST_PRESICION/dongdian/AIDD-3D-Pred/amber26/install/ambertools26/amber.sh


PROJ=/media/geno/3d3722c0-43ee-42ed-ae82-7a640071f8b1/ST_PRESICION/dongdian/AIDD-3D-Pred/6_epitopes_hla_pred

COCO=$PROJ/new_analysis/cocomaps


mkdir -p $COCO/pdb
mkdir -p $COCO/logs


declare -A systems

systems["E6wt"]="E6wt_10mer"
systems["E6apl"]="E6apl_10mer"
systems["E7wt1"]="E7wt1_9mer"
systems["E7apl1"]="E7apl1_9mer"
systems["E7wt2"]="E7wt2_10mer"
systems["E7apl2"]="E7apl2_10mer"


for name in "${!systems[@]}"
do

prefix=${systems[$name]}

echo "================================="
echo "Processing $name"
echo "================================="


cat > $COCO/scripts/${name}.in <<EOF

parm $PROJ/new_analysis/prmtops/${prefix}_com.prmtop

trajin $PROJ/analysis/clean_traj/${prefix}_clean.nc 4001 10000 10

autoimage

rms first :1-365@CA,C,N

trajout $COCO/pdb/${name}_cocomaps.pdb pdb onlyframes 1

run
quit

EOF


cpptraj -i $COCO/scripts/${name}.in \
> $COCO/logs/${name}.log 2>&1


done


echo ""
echo "================================="
echo "CoCoMaps2 PDB generation finished"
echo "================================="

ls -lh $COCO/pdb
