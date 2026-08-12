#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


input_file = Path(
    "../results/tables/epitope_conservation.tsv"
)


output_file = Path(
    "../results/tables/"
    "APL_anchor_conservation_Supplementary_Table.tsv"
)



epi = pd.read_csv(
    input_file,
    sep="\t"
)



mutation_info = [

[
"E6apl",
"T1Y",
"E6wt",
1,
"T",
"Y"
],

[
"E6apl",
"I2L",
"E6wt",
2,
"I",
"L"
],

[
"E7apl1",
"L1Y",
"E7wt1",
1,
"L",
"Y"
],

[
"E7apl2",
"T10V",
"E7wt2",
10,
"T",
"V"
],

[
"E7apl3",
"M2L",
"E7wt2",
2,
"M",
"L"
],

[
"E7apl3",
"T10V",
"E7wt2",
10,
"T",
"V"
]

]


results=[]


for apl,mutation,wtpep,pos,wt,aplres in mutation_info:


    row=epi[

        (epi.Peptide==wtpep)
        &
        (epi.Residue_position==pos)

    ]


    results.append({

        "APL":
        apl,

        "Mutation":
        mutation,

        "WT_residue":
        wt,

        "APL_residue":
        aplres,

        "Natural_major_residue":
        row.iloc[0]["Major_AA"],

        "WT_residue_conservation_percent":
        row.iloc[0]["Conservation_percent"]

    })



pd.DataFrame(results).to_csv(

    output_file,

    sep="\t",

    index=False

)


print(
"APL supplementary table generated"
)