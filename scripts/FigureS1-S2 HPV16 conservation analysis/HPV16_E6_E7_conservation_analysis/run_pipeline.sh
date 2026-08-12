#!/bin/bash

set -e

ROOT=$(pwd)

echo "Step1 download"
cd $ROOT/01_download
python download_hpv16_variants.py

test -f $ROOT/data/raw_nt/HPV16_E6_nt.fasta
test -f $ROOT/data/raw_nt/HPV16_E7_nt.fasta


echo "Step2 translation"
cd $ROOT/02_translation
python translate_E6_E7.py


test -f $ROOT/data/raw_aa/HPV16_E6_aa.fasta
test -f $ROOT/data/raw_aa/HPV16_E7_aa.fasta


echo "Step3 alignment"
cd $ROOT/03_alignment
bash run_mafft.sh


echo "Pipeline finished"