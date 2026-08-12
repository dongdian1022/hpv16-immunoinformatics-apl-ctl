#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Translate HPV16 E6/E7 nucleotide sequences into amino acid sequences.

Input:
    ../data/raw_nt/
        HPV16_E6_nt.fasta
        HPV16_E7_nt.fasta

Output:
    ../data/raw_aa/
        HPV16_E6_aa.fasta
        HPV16_E7_aa.fasta
"""


from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


# =========================
# Paths
# =========================

RAW_NT = Path("../data/raw_nt")
RAW_AA = Path("../data/raw_aa")


# create output directory automatically

RAW_AA.mkdir(
    parents=True,
    exist_ok=True
)


# =========================
# Translation function
# =========================

def translate_gene(gene):

    input_file = (
        RAW_NT /
        f"HPV16_{gene}_nt.fasta"
    )

    output_file = (
        RAW_AA /
        f"HPV16_{gene}_aa.fasta"
    )


    if not input_file.exists():

        raise FileNotFoundError(
            f"Missing input file: {input_file}"
        )


    translated_records = []


    for record in SeqIO.parse(
        input_file,
        "fasta"
    ):

        aa_seq = Seq(
            str(record.seq)
        ).translate(
            to_stop=True
        )


        record.seq = aa_seq

        record.description += (
            " translated_protein"
        )


        translated_records.append(
            record
        )


    SeqIO.write(
        translated_records,
        output_file,
        "fasta"
    )


    print(
        f"{gene}: "
        f"{len(translated_records)} sequences translated"
    )


# =========================
# Main
# =========================

if __name__ == "__main__":


    for gene in [
        "E6",
        "E7"
    ]:

        translate_gene(
            gene
        )


    print(
        "Translation completed."
    )