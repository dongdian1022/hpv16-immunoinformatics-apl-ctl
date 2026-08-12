#!/bin/bash
# ============================================================
# generate_all_chainAB.sh
# 用途：
#   1) 从 clean trajectory + com.prmtop 提取代表性PDB
#   2) 自动给PDB添加 chain ID
#      - HLA-A*02:01 = Chain A (1-365)
#      - peptide      = Chain B (376-384 或 376-385)
#
# 最终输出：
#   new_analysis/cocomaps/pdb/
#     E6wt_chainAB.pdb
#     E6apl_chainAB.pdb
#     E7wt1_chainAB.pdb
#     E7apl1_chainAB.pdb
#     E7wt2_chainAB.pdb
#     E7apl2_chainAB.pdb
# ============================================================

set -e

source /ST_PRESICION/dongdian/AIDD-3D-Pred/amber26/install/ambertools26/amber.sh

PROJ=/media/geno/3d3722c0-43ee-42ed-ae82-7a640071f8b1/ST_PRESICION/dongdian/AIDD-3D-Pred/6_epitopes_hla_pred
COCO=$PROJ/new_analysis/cocomaps
PRM=$PROJ/new_analysis/prmtops
TRJ=$PROJ/analysis/clean_traj

mkdir -p "$COCO/scripts" "$COCO/pdb" "$COCO/logs"

echo "==================================================="
echo "🚀 批量生成 6 个 CoCoMaps2 上传用 chainAB PDB"
echo "项目目录: $PROJ"
echo "输出目录: $COCO/pdb"
echo "开始时间: $(date)"
echo "==================================================="

# 体系列表： 名称 前缀 peptide终止残基号
SYSTEMS=(
"E6wt   E6wt_10mer   385"
"E6apl  E6apl_10mer  385"
"E7wt1  E7wt1_9mer   384"
"E7apl1 E7apl1_9mer  384"
"E7wt2  E7wt2_10mer  385"
"E7apl2 E7apl2_10mer 385"
)

for item in "${SYSTEMS[@]}"; do
    name=$(echo $item | awk '{print $1}')
    prefix=$(echo $item | awk '{print $2}')
    pep_end=$(echo $item | awk '{print $3}')

    raw_pdb="$COCO/pdb/${name}_raw.pdb"
    out_pdb="$COCO/pdb/${name}_chainAB.pdb"
    in_file="$COCO/scripts/${name}_export.in"
    log_file="$COCO/logs/${name}.log"

    echo ""
    echo "---------------------------------------------------"
    echo ">> 处理体系: $name"
    echo "   prefix   = $prefix"
    echo "   pep_end  = $pep_end"
    echo "---------------------------------------------------"

    # ---------- Step 1: cpptraj 导出 raw pdb ----------
    cat > "$in_file" << EOF
parm $PRM/${prefix}_com.prmtop

trajin $TRJ/${prefix}_clean.nc 4001 10000 10

autoimage

rms first :1-365@CA,C,N

trajout $raw_pdb pdb onlyframes 1

run
quit
EOF

    cpptraj -i "$in_file" > "$log_file" 2>&1

    if [ ! -s "$raw_pdb" ]; then
        echo "❌ raw pdb 生成失败: $raw_pdb"
        echo "   请查看日志: $log_file"
        exit 1
    fi

    # ---------- Step 2: 添加 chain ID ----------
    python3 - << PY
input_pdb = r"$raw_pdb"
output_pdb = r"$out_pdb"
pep_end = $pep_end

with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM", "HETATM")):
            try:
                resnum = int(line[22:26])
            except:
                fout.write(line)
                continue

            if 1 <= resnum <= 365:
                chain = "A"
            elif 376 <= resnum <= pep_end:
                chain = "B"
            else:
                chain = " "

            line = line[:21] + chain + line[22:]

        fout.write(line)

print("Generated:", output_pdb)
PY

    if [ ! -s "$out_pdb" ]; then
        echo "❌ chainAB pdb 生成失败: $out_pdb"
        exit 1
    fi

    echo "✅ 完成: $out_pdb"
done

echo ""
echo "==================================================="
echo "🎉 全部完成！最终可上传 CoCoMaps2 的文件如下："
ls -lh "$COCO/pdb"/*_chainAB.pdb
echo "结束时间: $(date)"
echo "==================================================="
