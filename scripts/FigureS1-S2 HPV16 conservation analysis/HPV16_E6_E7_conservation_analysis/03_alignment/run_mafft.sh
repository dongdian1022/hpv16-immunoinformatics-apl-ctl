#!/bin/bash

mafft --localpair --maxiterate 1000 ../data/raw_aa/HPV16_E6_aa.fasta > ../data/aligned/HPV16_E6_AA_aligned.fasta

mafft --localpair --maxiterate 1000 ../data/raw_aa/HPV16_E7_aa.fasta > ../data/aligned/HPV16_E7_AA_aligned.fasta