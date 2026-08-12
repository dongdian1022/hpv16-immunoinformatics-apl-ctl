#!/usr/bin/env python3

input_pdb = "pdb/E7apl1_raw.pdb"
output_pdb = "pdb/E7apl1_chainAB.pdb"


with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:

    for line in fin:

        if line.startswith(("ATOM", "HETATM")):

            try:
                resnum = int(line[22:26])
            except:
                fout.write(line)
                continue

            # HLA-A*02:01
            if 1 <= resnum <= 365:
                chain = "A"

            # peptide YLMGTLGIV
            elif 376 <= resnum <= 384:
                chain = "B"

            else:
                chain = " "

            # PDB chain位置：第22列(index 21)
            line = line[:21] + chain + line[22:]

        fout.write(line)


print("Finished:")
print(output_pdb)
