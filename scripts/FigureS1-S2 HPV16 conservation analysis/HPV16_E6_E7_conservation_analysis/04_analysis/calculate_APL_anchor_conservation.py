#!/usr/bin/env python3


from pathlib import Path
import pandas as pd



OUT = Path(
    "../results/tables"
)


OUT.mkdir(
    parents=True,
    exist_ok=True
)



data=[

[
"E6apl",
"T1Y",
"TIHDIILECV",
"YLHDIILECV"
],

[
"E6apl",
"I2L",
"TIHDIILECV",
"YLHDIILECV"
],

[
"E7apl1",
"L1Y",
"LLMGTLGIV",
"YLMGTLGIV"
],

[
"E7apl2",
"T10V",
"YMLDLQPETT",
"YMLDLQPETV"
],

[
"E7apl3",
"M2L",
"YMLDLQPETT",
"YLLDLQPETV"
],

[
"E7apl3",
"T10V",
"YMLDLQPETT",
"YLLDLQPETV"

]

]


df=pd.DataFrame(

    data,

    columns=[
        "APL",
        "Mutation",
        "WT_sequence",
        "APL_sequence"
    ]

)



df["Conservation_status"]="pending"


df.to_csv(

    OUT /
    "APL_anchor_conservation.tsv",

    sep="\t",

    index=False

)


print(
    "APL table generated"
)