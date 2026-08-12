#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


table_dir=Path(
    "../results/tables"
)


files=list(
    table_dir.glob("*.tsv")
)


summary=[]


for f in files:

    df=pd.read_csv(
        f,
        sep="\t"
    )


    summary.append({

        "File":
        f.name,

        "Rows":
        len(df)

    })


pd.DataFrame(summary).to_csv(

    table_dir /
    "analysis_summary.tsv",

    sep="\t",

    index=False

)


print(
    "Summary generated"
)