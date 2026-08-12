#!/bin/bash
# =============================================================================
# 05_extract_mmgbsa.sh
# Parse MMPBSA.py output files into one summary table + CSV.
# Columns: System, TOTAL, VDWAALS, EEL, EGB, ESURF (all kcal/mol, mean +/- SD).
# 将 MMPBSA 结果汇总为三线表 / CSV。
#
# Input : analysis/mmpbsa/data/*_FINAL_energy.dat
# Output: analysis/mmpbsa/data/MMGBSA_Summary.csv (+ printed table)
# Usage : bash 05_extract_mmgbsa.sh
# =============================================================================
set -euo pipefail

data_dir="analysis/mmpbsa/data"
csv_file="${data_dir}/MMGBSA_Summary.csv"
echo "System,TOTAL(kcal/mol),VDWAALS,EEL,EGB,ESURF" > "$csv_file"

printf "%-16s | %-14s | %-14s | %-14s | %-14s | %-14s\n" \
       "System" "TOTAL" "VDWAALS" "EEL" "EGB" "ESURF"
printf -- '-%.0s' {1..100}; echo

shopt -s nullglob
for file in "$data_dir"/*_FINAL_energy.dat; do
    name=$(basename "$file" _FINAL_energy.dat)
    awk -v sys="$name" -v csv="$csv_file" '
    BEGIN {
        vdw="N/A"; eel="N/A"; egb="N/A"; esurf="N/A"; tot="N/A";
        vdw_c="N/A"; eel_c="N/A"; egb_c="N/A"; esurf_c="N/A"; tot_c="N/A"
    }
    /Differences/ {flag=1}
    flag && $1=="VDWAALS" {vdw=sprintf("%.2f ± %.2f",$2,$3); vdw_c=sprintf("%.2f±%.2f",$2,$3)}
    flag && $1=="EEL"     {eel=sprintf("%.2f ± %.2f",$2,$3); eel_c=sprintf("%.2f±%.2f",$2,$3)}
    flag && $1=="EGB"     {egb=sprintf("%.2f ± %.2f",$2,$3); egb_c=sprintf("%.2f±%.2f",$2,$3)}
    flag && $1=="ESURF"   {esurf=sprintf("%.2f ± %.2f",$2,$3); esurf_c=sprintf("%.2f±%.2f",$2,$3)}
    # total binding free energy (covers AmberTools naming variants)
    flag && $1=="DELTA" && $2=="G" && $3=="bind" {tot=sprintf("%.2f ± %.2f",$4,$5); tot_c=sprintf("%.2f±%.2f",$4,$5)}
    flag && $1=="DELTA" && $2=="TOTAL"           {tot=sprintf("%.2f ± %.2f",$3,$4); tot_c=sprintf("%.2f±%.2f",$3,$4)}
    flag && $1=="TOTAL" && $2!="Energy" && tot=="N/A" {tot=sprintf("%.2f ± %.2f",$2,$3); tot_c=sprintf("%.2f±%.2f",$2,$3)}
    END {
        printf "%-16s | %-14s | %-14s | %-14s | %-14s | %-14s\n", sys, tot, vdw, eel, egb, esurf
        print sys","tot_c","vdw_c","eel_c","egb_c","esurf_c >> csv
    }' "$file"
done

printf -- '-%.0s' {1..100}; echo
echo "CSV written to: $csv_file"
