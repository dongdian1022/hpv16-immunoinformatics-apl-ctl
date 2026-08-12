#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01_extract_and_prepare_cocomaps2.py

Purpose
-------
Extract CoCoMaps2 result zip archives and build all summary tables required for
Figure 8 and Supplementary CoCoMaps2 figures.

Input
-----
raw_zip/*.zip

Output
------
summary/
    cocomaps_master_summary.csv
    cocomaps_interaction_matrix.csv
    cocomaps_peptide_residue_contacts.csv
    cocomaps_anchor_contact.csv
    cocomaps_hla_hotspots.csv
    cocomaps_pairwise_delta.csv
"""

import os
import re
import zipfile
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd


# =========================
# Paths
# =========================
ROOT = "/ST_PRESICION/dongdian/AIDD-3D-Pred/cocomaps2/"
ZIP_DIR = os.path.join(ROOT, "raw_zip")
OUT_DIR = os.path.join(ROOT, "summary")
os.makedirs(OUT_DIR, exist_ok=True)


# =========================
# System metadata
# =========================
SYSTEM_META = OrderedDict({
    "E6wt":   {"pair": "E6",    "type": "WT",  "seq": "TIHDIILECV", "length": 10},
    "E6apl":  {"pair": "E6",    "type": "APL", "seq": "YLHDIILECV", "length": 10},
    "E7wt1":  {"pair": "E7-1",  "type": "WT",  "seq": "LLMGTLGIV",  "length": 9},
    "E7apl1": {"pair": "E7-1",  "type": "APL", "seq": "YLMGTLGIV",  "length": 9},
    "E7wt2":  {"pair": "E7-2",  "type": "WT",  "seq": "YMLDLQPETT", "length": 10},
    "E7apl2": {"pair": "E7-2",  "type": "APL", "seq": "YMLDLQPETV", "length": 10},
})

SYSTEM_ORDER = ["E6wt", "E6apl", "E7wt1", "E7apl1", "E7wt2", "E7apl2"]

PAIRINGS = [
    ("E6wt", "E6apl", "E6"),
    ("E7wt1", "E7apl1", "E7-1"),
    ("E7wt2", "E7apl2", "E7-2"),
]

INTERACTION_FILE_MAP = {
    "H-bond": "H_bond",
    "Apolar_vdw": "Apolar_vdw",
    "Polar_vdw": "Polar_vdw",
    "C-H_ON": "CH_O_N",
    "C-H_pi": "CH_pi",
    "pi-pi": "pi_pi",
    "Amino_pi": "Amino_pi",
    "Anion_pi": "Anion_pi",
    "Cation_pi": "Cation_pi",
    "Proximal": "Proximal",
    "Clash": "Clash",
    "Salt_bridge": "Salt_bridge",
    "Water_Mediated": "Water_Mediated",
    "Metal_Mediated": "Metal_Mediated",
    "Halogen_bond": "Halogen_bond",
    "SS_bond": "SS_bond",
    "Lone_pair_pi": "Lone_pair_pi",
    "N-S-O-H_pi": "NSOH_pi"
}

SELECTED_INTERACTIONS = [
    "H_bond",
    "CH_O_N",
    "Apolar_vdw",
    "Polar_vdw",
    "CH_pi",
    "pi_pi",
    "Amino_pi",
    "Anion_pi",
    "Cation_pi",
    "Proximal",
    "Clash"
]


# =========================
# Helper functions
# =========================
def detect_system(text):
    text_low = text.lower()
    for sys_name in SYSTEM_META.keys():
        if sys_name.lower() in text_low:
            return sys_name
    return None


def safe_read_csv_from_zip(zf, contains_text):
    for name in zf.namelist():
        if contains_text in name and name.lower().endswith(".csv"):
            try:
                with zf.open(name) as f:
                    return pd.read_csv(f)
            except Exception:
                return None
    return None


def safe_count_rows_from_zip(zf, internal_name):
    try:
        with zf.open(internal_name) as f:
            df = pd.read_csv(f)
        return df.shape[0]
    except Exception:
        return 0


def parse_area_pair(value):
    """
    Parse strings like:
    '1785.1 / 892.55'
    """
    try:
        parts = [x.strip() for x in str(value).split("/")]
        if len(parts) == 2:
            return float(parts[0]), float(parts[1])
    except Exception:
        pass
    return np.nan, np.nan


def detect_interaction_types(type_str):
    """
    Parse the 'Type of Interactions' text in final_file.csv
    """
    s = str(type_str)
    found = set()

    if "H-bond" in s:
        found.add("H_bond")
    if "CH-O/N bond" in s:
        found.add("CH_O_N")
    if "Apolar vdW contact" in s:
        found.add("Apolar_vdw")
    if "Polar vdW contact" in s:
        found.add("Polar_vdw")
    if "CH-pi" in s or "CH-pi interaction" in s or "C-H pi" in s:
        found.add("CH_pi")
    if "pi-pi" in s:
        found.add("pi_pi")
    if "Amino-pi" in s or "Amino_pi" in s:
        found.add("Amino_pi")
    if "Anion-pi" in s or "Anion_pi" in s:
        found.add("Anion_pi")
    if "Cation-pi" in s or "Cation_pi" in s:
        found.add("Cation_pi")
    if "Proximal" in s:
        found.add("Proximal")
    if "Clash" in s:
        found.add("Clash")

    return found


def reorder_systems(df, col="System"):
    df[col] = pd.Categorical(df[col], categories=SYSTEM_ORDER, ordered=True)
    return df.sort_values(col).reset_index(drop=True)


# =========================
# Main extraction
# =========================
zip_files = sorted([os.path.join(ZIP_DIR, x) for x in os.listdir(ZIP_DIR) if x.endswith(".zip")])

if len(zip_files) == 0:
    raise FileNotFoundError(f"❌ No zip files found in: {ZIP_DIR}")

master_rows = []
peptide_rows = []
hotspot_rows = []

for zip_path in zip_files:
    with zipfile.ZipFile(zip_path, "r") as zf:
        joined_names = " ".join(zf.namelist()) + " " + os.path.basename(zip_path)
        system = detect_system(joined_names)

        if system is None:
            print(f"⚠️ Cannot detect system for: {os.path.basename(zip_path)} ; skipped")
            continue

        meta = SYSTEM_META[system]

        # ---------- count interaction-category files ----------
        interaction_counts = defaultdict(int)

        for name in zf.namelist():
            if not name.lower().endswith(".csv"):
                continue

            base = os.path.basename(name)

            if (
                "final_file.csv" in base
                or "small_summary.csv" in base
                or "summary_table.csv" in base
                or "ASA_table_chain1.csv" in base
                or "ASA_table_chain2.csv" in base
                or "Rsa_stats.csv" in base
            ):
                continue

            m = re.search(r"_A_B_(.+)\.csv$", base)
            if not m:
                continue

            raw_name = m.group(1)
            label = INTERACTION_FILE_MAP.get(raw_name, raw_name)

            interaction_counts[label] = safe_count_rows_from_zip(zf, name)

        # ---------- small summary ----------
        small_summary = safe_read_csv_from_zip(zf, "small_summary.csv")
        total_interacting_residues = np.nan
        chain1_interacting_residues = np.nan
        chain2_interacting_residues = np.nan

        if small_summary is not None and small_summary.shape[1] >= 3:
            prop_col = small_summary.columns[1]
            val_col = small_summary.columns[2]
            mapping = dict(zip(small_summary[prop_col], small_summary[val_col]))

            total_interacting_residues = mapping.get("Total Number of Interacting Residues", np.nan)
            chain1_interacting_residues = mapping.get("Number of Interacting Residues in Chain 1", np.nan)
            chain2_interacting_residues = mapping.get("Number of Interacting Residues in Chain 2", np.nan)

        # ---------- RSA / interface area ----------
        rsa_stats = safe_read_csv_from_zip(zf, "Rsa_stats.csv")
        buried_area = np.nan
        interface_area = np.nan
        polar_buried_area = np.nan
        polar_interface_area = np.nan
        nonpolar_buried_area = np.nan
        nonpolar_interface_area = np.nan

        if rsa_stats is not None and rsa_stats.shape[1] >= 3:
            prop_col = rsa_stats.columns[1]
            val_col = rsa_stats.columns[2]
            rsa_map = dict(zip(rsa_stats[prop_col], rsa_stats[val_col]))

            buried_area, interface_area = parse_area_pair(
                rsa_map.get("Buried area upon the complex formation / Interface area (Å²)", np.nan)
            )
            polar_buried_area, polar_interface_area = parse_area_pair(
                rsa_map.get("POLAR Buried area upon the complex formation / Interface area (Å²)", np.nan)
            )
            nonpolar_buried_area, nonpolar_interface_area = parse_area_pair(
                rsa_map.get("NON POLAR Buried area upon the complex formation / Interface area (Å²)", np.nan)
            )

        # ---------- ASA chain-specific ----------
        asa1 = safe_read_csv_from_zip(zf, "ASA_table_chain1.csv")
        asa2 = safe_read_csv_from_zip(zf, "ASA_table_chain2.csv")

        chain1_buried_sum = asa1["Buried ASA (Interface)"].sum() if asa1 is not None and "Buried ASA (Interface)" in asa1.columns else np.nan
        chain2_buried_sum = asa2["Buried ASA (Interface)"].sum() if asa2 is not None and "Buried ASA (Interface)" in asa2.columns else np.nan

        # ---------- final_file.csv ----------
        final_df = safe_read_csv_from_zip(zf, "final_file.csv")
        residue_pairs = np.nan

        pep_contact_counter = defaultdict(int)
        pep_type_counter = defaultdict(lambda: defaultdict(int))

        hla_contact_counter = defaultdict(int)
        hla_type_counter = defaultdict(lambda: defaultdict(int))

        if final_df is not None:
            final_df.columns = [str(x).strip() for x in final_df.columns]
            residue_pairs = final_df.shape[0]

            for _, row in final_df.iterrows():
                try:
                    chain1 = str(row["Chain 1"]).strip()
                    chain2 = str(row["Chain 2"]).strip()

                    res1_name = str(row["Res. Name 1"]).strip()
                    res1_num  = int(row["Res. Number 1"])

                    res2_name = str(row["Res. Name 2"]).strip()
                    res2_num  = int(row["Res. Number 2"])

                    types = detect_interaction_types(row["Type of Interactions"])
                except Exception:
                    continue

                if chain1 == "A" and chain2 == "B":
                    pep_pos = res2_num - 375
                    pep_key = (pep_pos, res2_num, res2_name)
                    hla_key = (res1_num, res1_name)

                    pep_contact_counter[pep_key] += 1
                    hla_contact_counter[hla_key] += 1

                    for t in types:
                        pep_type_counter[pep_key][t] += 1
                        hla_type_counter[hla_key][t] += 1

        # ---------- peptide residue contribution ----------
        for (pep_pos, pep_resnum, pep_resname), cnt in sorted(pep_contact_counter.items(), key=lambda x: x[0][0]):
            out = {
                "System": system,
                "Pair": meta["pair"],
                "Type": meta["type"],
                "Sequence": meta["seq"],
                "PeptideLength": meta["length"],
                "PeptidePosition": pep_pos,
                "PeptideResidueNumber": pep_resnum,
                "PeptideResidueName": pep_resname,
                "ContactPairs": cnt
            }
            for t in SELECTED_INTERACTIONS:
                out[t] = pep_type_counter[(pep_pos, pep_resnum, pep_resname)].get(t, 0)
            peptide_rows.append(out)

        # ---------- HLA hotspot ----------
        for (resnum, resname), cnt in sorted(hla_contact_counter.items(), key=lambda x: x[0][0]):
            out = {
                "System": system,
                "Pair": meta["pair"],
                "Type": meta["type"],
                "HLAResidueNumber": resnum,
                "HLAResidueName": resname,
                "HLAResidue": f"{resname.upper()}{resnum}",
                "ContactPairs": cnt
            }
            for t in SELECTED_INTERACTIONS:
                out[t] = hla_type_counter[(resnum, resname)].get(t, 0)
            hotspot_rows.append(out)

        favorable_contacts = sum([
            interaction_counts.get("H_bond", 0),
            interaction_counts.get("CH_O_N", 0),
            interaction_counts.get("Apolar_vdw", 0),
            interaction_counts.get("Polar_vdw", 0),
            interaction_counts.get("CH_pi", 0),
            interaction_counts.get("pi_pi", 0),
            interaction_counts.get("Amino_pi", 0),
            interaction_counts.get("Anion_pi", 0),
            interaction_counts.get("Cation_pi", 0),
            interaction_counts.get("Proximal", 0)
        ])

        row = {
            "System": system,
            "Pair": meta["pair"],
            "Type": meta["type"],
            "Sequence": meta["seq"],
            "PeptideLength": meta["length"],
            "ResiduePairs": residue_pairs,
            "TotalInteractingResidues": total_interacting_residues,
            "Chain1InteractingResidues": chain1_interacting_residues,
            "Chain2InteractingResidues": chain2_interacting_residues,
            "BuriedArea_A2": buried_area,
            "InterfaceArea_A2": interface_area,
            "PolarBuriedArea_A2": polar_buried_area,
            "PolarInterfaceArea_A2": polar_interface_area,
            "NonPolarBuriedArea_A2": nonpolar_buried_area,
            "NonPolarInterfaceArea_A2": nonpolar_interface_area,
            "Chain1BuriedASA_Sum": chain1_buried_sum,
            "Chain2BuriedASA_Sum": chain2_buried_sum,
            "FavorableContacts": favorable_contacts,
        }

        for label in SELECTED_INTERACTIONS:
            row[label] = interaction_counts.get(label, 0)

        master_rows.append(row)


# =========================
# Build output tables
# =========================
master_df = pd.DataFrame(master_rows)
master_df = reorder_systems(master_df, "System")

interaction_matrix_df = master_df[[
    "System", "Pair", "Type",
    "H_bond", "CH_O_N", "Apolar_vdw", "Polar_vdw",
    "CH_pi", "pi_pi", "Amino_pi", "Anion_pi",
    "Cation_pi", "Proximal", "Clash", "FavorableContacts"
]].copy()

peptide_df = pd.DataFrame(peptide_rows)
peptide_df = reorder_systems(peptide_df, "System")

hotspot_df = pd.DataFrame(hotspot_rows)
hotspot_df = reorder_systems(hotspot_df, "System")

# Anchor contact table
anchor_rows = []
for sys_name, meta in SYSTEM_META.items():
    sub = peptide_df.loc[peptide_df["System"] == sys_name].copy()
    if sub.empty:
        continue

    cterm_pos = int(sub["PeptidePosition"].max())

    for label, pos in [("P1", 1), ("P2", 2), ("Cterm", cterm_pos)]:
        tmp = sub.loc[sub["PeptidePosition"] == pos]
        if len(tmp) == 0:
            cnt = 0
            resname = ""
        else:
            cnt = int(tmp["ContactPairs"].iloc[0])
            resname = str(tmp["PeptideResidueName"].iloc[0])

        anchor_rows.append({
            "System": sys_name,
            "Pair": meta["pair"],
            "Type": meta["type"],
            "AnchorPosition": label,
            "PeptidePosition": pos,
            "ResidueName": resname,
            "ContactPairs": cnt
        })

anchor_df = pd.DataFrame(anchor_rows)
anchor_df = reorder_systems(anchor_df, "System")

# Pairwise delta
delta_rows = []
compare_cols = [
    "ResiduePairs", "TotalInteractingResidues",
    "BuriedArea_A2", "InterfaceArea_A2",
    "H_bond", "CH_O_N", "Apolar_vdw", "Polar_vdw",
    "CH_pi", "pi_pi", "Amino_pi", "Anion_pi",
    "Cation_pi", "Proximal", "Clash", "FavorableContacts"
]

for wt_sys, apl_sys, pair_name in PAIRINGS:
    wt = master_df.loc[master_df["System"] == wt_sys].iloc[0]
    apl = master_df.loc[master_df["System"] == apl_sys].iloc[0]

    out = {
        "Pair": pair_name,
        "WT_System": wt_sys,
        "APL_System": apl_sys,
        "WT_Sequence": wt["Sequence"],
        "APL_Sequence": apl["Sequence"]
    }

    for c in compare_cols:
        out[f"{c}_WT"] = wt[c]
        out[f"{c}_APL"] = apl[c]
        out[f"{c}_Delta_APL_minus_WT"] = apl[c] - wt[c]

    delta_rows.append(out)

delta_df = pd.DataFrame(delta_rows)

# Save
master_df.to_csv(os.path.join(OUT_DIR, "cocomaps_master_summary.csv"), index=False)
interaction_matrix_df.to_csv(os.path.join(OUT_DIR, "cocomaps_interaction_matrix.csv"), index=False)
peptide_df.to_csv(os.path.join(OUT_DIR, "cocomaps_peptide_residue_contacts.csv"), index=False)
anchor_df.to_csv(os.path.join(OUT_DIR, "cocomaps_anchor_contact.csv"), index=False)
hotspot_df.to_csv(os.path.join(OUT_DIR, "cocomaps_hla_hotspots.csv"), index=False)
delta_df.to_csv(os.path.join(OUT_DIR, "cocomaps_pairwise_delta.csv"), index=False)

print("✅ Finished.")
print(os.path.join(OUT_DIR, "cocomaps_master_summary.csv"))
print(os.path.join(OUT_DIR, "cocomaps_interaction_matrix.csv"))
print(os.path.join(OUT_DIR, "cocomaps_peptide_residue_contacts.csv"))
print(os.path.join(OUT_DIR, "cocomaps_anchor_contact.csv"))
print(os.path.join(OUT_DIR, "cocomaps_hla_hotspots.csv"))
print(os.path.join(OUT_DIR, "cocomaps_pairwise_delta.csv"))