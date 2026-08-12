#!/usr/bin/env python3

from pathlib import Path
from Bio import AlignIO
import pandas as pd


ALIGN_DIR = Path("../data/aligned")
OUT_DIR = Path("../results/tables")

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)



epitopes = {

"E6wt":
("E6",29,"TIHDIILECV"),

"E6apl":
("E6",29,"YLHDIILECV"),

"E7wt1":
("E7",82,"LLMGTLGIV"),

"E7apl1":
("E7",82,"YLMGTLGIV"),

"E7wt2":
("E7",11,"YMLDLQPETT"),

"E7apl2":
("E7",11,"YMLDLQPETV"),

"E7apl3":
("E7",11,"YLLDLQPETV")

}



rows=[]



for peptide,(gene,start,sequence) in epitopes.items():


    alignment = AlignIO.read(

        ALIGN_DIR /
        f"HPV16_{gene}_AA_aligned.fasta",

        "fasta"

    )


    sequences = [
        str(record.seq)
        for record in alignment
    ]


    for i,aa in enumerate(
        sequence
    ):


        # reference coordinate
        index = start+i-1


        observed=[]


        for seq in sequences:

            if index < len(seq):

                if seq[index]!="-" :

                    observed.append(
                        seq[index]
                    )


        if len(observed)==0:
            continue


        major=max(
            set(observed),
            key=observed.count
        )


        conservation=(
            observed.count(major)
            /
            len(observed)
            *
            100
        )


        rows.append({

            "Peptide":
            peptide,

            "Residue_position":
            i+1,

            "Reference_AA":
            aa,

            "Major_AA":
            major,

            "Conservation_percent":
            round(
                conservation,
                2
            )

        })



df=pd.DataFrame(rows)



df.to_csv(

    OUT_DIR /
    "epitope_conservation.tsv",

    sep="\t",

    index=False

)


print(
    "Epitope conservation completed"
)