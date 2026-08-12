#!/usr/bin/env python3

from pathlib import Path
import time
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = "your_email@example.com"

ACC = Path("../config/hpv16_accessions.tsv")
OUT = Path("../data/raw_nt")
META = Path("../results/tables")
LOG = Path("../results/logs")

OUT.mkdir(parents=True, exist_ok=True)
META.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)


def fetch_record(accession, retry=3):
    for i in range(retry):
        try:
            h = Entrez.efetch(
                db="nuccore",
                id=accession,
                rettype="gb",
                retmode="text"
            )
            r = SeqIO.read(h, "genbank")
            h.close()
            return r
        except Exception:
            time.sleep(2)
    return None


def extract_gene(record, gene):
    gene = gene.upper()

    for feature in record.features:

        if feature.type != "CDS":
            continue

        text=[]

        for key in [
            "gene",
            "product",
            "note",
            "protein_id"
        ]:
            text.extend(
                feature.qualifiers.get(key, [])
            )

        annotation=" ".join(text).upper()

        if gene in annotation:
            return feature.extract(record.seq)

    return None


def main():

    table=pd.read_csv(
        ACC,
        sep="\t"
    )

    e6=[]
    e7=[]
    qc=[]

    for _,row in table.iterrows():

        record=fetch_record(
            row.Accession
        )

        if record is None:
            qc.append({
                "Lineage":row.Lineage,
                "Accession":row.Accession,
                "Status":"download_failed"
            })
            continue


        seq6=extract_gene(record,"E6")
        seq7=extract_gene(record,"E7")


        if seq6:
            e6.append(
                SeqRecord(
                    seq6,
                    id=row.Lineage,
                    description=row.Accession
                )
            )

        if seq7:
            e7.append(
                SeqRecord(
                    seq7,
                    id=row.Lineage,
                    description=row.Accession
                )
            )


        qc.append({
            "Lineage":row.Lineage,
            "Accession":row.Accession,
            "E6_found":bool(seq6),
            "E7_found":bool(seq7)
        })


    SeqIO.write(
        e6,
        OUT/"HPV16_E6_nt.fasta",
        "fasta"
    )

    SeqIO.write(
        e7,
        OUT/"HPV16_E7_nt.fasta",
        "fasta"
    )


    pd.DataFrame(qc).to_csv(
        META/"HPV16_variant_QC.tsv",
        sep="\t",
        index=False
    )


if __name__=="__main__":
    main()