#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from Bio import AlignIO
import pandas as pd


ALIGN_DIR = Path("../data/aligned")
OUT_DIR = Path("../results/tables")

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)


def calculate_conservation(
        alignment_file,
        gene
):

    alignment = AlignIO.read(
        alignment_file,
        "fasta"
    )

    sequences = [
        str(record.seq)
        for record in alignment
    ]


    results = []


    for pos in range(
        alignment.get_alignment_length()
    ):

        column = [
            seq[pos]
            for seq in sequences
        ]


        column = [
            aa
            for aa in column
            if aa != "-"
        ]


        if len(column)==0:
            continue


        major_aa = max(
            set(column),
            key=column.count
        )


        conservation = (
            column.count(major_aa)
            /
            len(column)
            *
            100
        )


        results.append({

            "Gene":
            gene,

            "Alignment_position":
            pos+1,

            "Major_AA":
            major_aa,

            "Conservation_percent":
            round(
                conservation,
                2
            ),

            "N_sequences":
            len(column)

        })


    return pd.DataFrame(results)



for gene in [
    "E6",
    "E7"
]:

    df = calculate_conservation(

        ALIGN_DIR /
        f"HPV16_{gene}_AA_aligned.fasta",

        gene

    )


    df.to_csv(

        OUT_DIR /
        f"{gene}_full_length_conservation.tsv",

        sep="\t",

        index=False

    )


    print(
        gene,
        "completed"
    )