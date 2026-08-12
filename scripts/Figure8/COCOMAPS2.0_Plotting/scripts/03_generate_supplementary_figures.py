#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_generate_CoCoMaps2_supplementary.py

Generate CoCoMaps2 supplementary figures.

Panels
------
A. Interaction-type heatmap
B. Peptide residue contribution heatmap
C. HLA hotspot heatmap
D. WT/APL delta heatmap

Outputs
-------
figures/Supplementary/
    Supplementary_CoCoMaps2_all.png
    Supplementary_CoCoMaps2_all.pdf
    Supplementary_CoCoMaps2_all.svg
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams


# =========================
# Paths
# =========================
ROOT = "/ST_PRESICION/dongdian/AIDD-3D-Pred/cocomaps2"
SUMMARY_DIR = os.path.join(ROOT, "summary")
OUT_DIR = os.path.join(ROOT, "figures", "Supplementary")
os.makedirs(OUT_DIR, exist_ok=True)


# =========================
# Style
# =========================
rcParams["font.family"] = "Arial"
rcParams["font.size"] = 9
rcParams["axes.linewidth"] = 1.0
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42

SYSTEM_ORDER = ["E6wt", "E6apl", "E7wt1", "E7apl1", "E7wt2", "E7apl2"]


# =========================
# Load
# =========================
master = pd.read_csv(os.path.join(SUMMARY_DIR, "cocomaps_master_summary.csv"))
peptide = pd.read_csv(os.path.join(SUMMARY_DIR, "cocomaps_peptide_residue_contacts.csv"))
hotspot = pd.read_csv(os.path.join(SUMMARY_DIR, "cocomaps_hla_hotspots.csv"))
delta = pd.read_csv(os.path.join(SUMMARY_DIR, "cocomaps_pairwise_delta.csv"))

master["System"] = pd.Categorical(master["System"], categories=SYSTEM_ORDER, ordered=True)
master = master.sort_values("System").reset_index(drop=True)

peptide["System"] = pd.Categorical(peptide["System"], categories=SYSTEM_ORDER, ordered=True)
peptide = peptide.sort_values(["System", "PeptidePosition"]).reset_index(drop=True)

hotspot["System"] = pd.Categorical(hotspot["System"], categories=SYSTEM_ORDER, ordered=True)
hotspot = hotspot.sort_values(["System", "HLAResidueNumber"]).reset_index(drop=True)


# =========================
# Helpers
# =========================
def panel_label(ax, label):
    ax.text(-0.12, 1.08, label, transform=ax.transAxes,
            fontsize=13, fontweight="bold", va="top", ha="left")


def draw_box_spines(ax):
    for side in ["top", "right", "bottom", "left"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(1.0)


def draw_heatmap(ax, data, xlabels, ylabels, title, cmap):
    im = ax.imshow(data, aspect="auto", cmap=cmap)

    ax.set_xticks(np.arange(len(xlabels)))
    ax.set_xticklabels(xlabels, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(ylabels)))
    ax.set_yticklabels(ylabels)
    ax.set_title(title)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=7)

    draw_box_spines(ax)
    return im


# =========================
# Panel A: Interaction-type heatmap
# =========================
interaction_cols = [
    "H_bond", "CH_O_N", "Apolar_vdw", "Polar_vdw",
    "CH_pi", "pi_pi", "Amino_pi", "Anion_pi",
    "Cation_pi", "Proximal", "Clash"
]

inter_mat = master[interaction_cols].to_numpy(dtype=float)


# =========================
# Panel B: Peptide residue contact heatmap
# =========================
pep_wide = peptide.pivot_table(
    index="System",
    columns="PeptidePosition",
    values="ContactPairs",
    aggfunc="sum"
).reindex(SYSTEM_ORDER)

pep_wide = pep_wide.reindex(columns=list(range(1, 11)))
pep_mat = pep_wide.to_numpy(dtype=float)
pep_xlabels = [f"P{x}" for x in pep_wide.columns]


# =========================
# Panel C: HLA hotspot heatmap
# top 15 hotspots
# =========================
top_hotspots = (
    hotspot.groupby("HLAResidue", as_index=False)["ContactPairs"]
    .sum()
    .sort_values("ContactPairs", ascending=False)
    .head(15)["HLAResidue"]
    .tolist()
)

hot_wide = hotspot.pivot_table(
    index="System",
    columns="HLAResidue",
    values="ContactPairs",
    aggfunc="sum"
).reindex(SYSTEM_ORDER)

hot_wide = hot_wide.reindex(columns=top_hotspots)
hot_mat = hot_wide.to_numpy(dtype=float)


# =========================
# Panel D: pairwise delta heatmap
# =========================
delta_metrics = [
    "InterfaceArea_A2_Delta_APL_minus_WT",
    "H_bond_Delta_APL_minus_WT",
    "CH_O_N_Delta_APL_minus_WT",
    "Apolar_vdw_Delta_APL_minus_WT",
    "Polar_vdw_Delta_APL_minus_WT",
    "FavorableContacts_Delta_APL_minus_WT"
]

delta_mat = delta[delta_metrics].to_numpy(dtype=float)
delta_ylabels = delta["Pair"].tolist()
delta_xlabels = [
    "ΔInterface area",
    "ΔH-bond",
    "ΔCH-O/N",
    "ΔApolar vdW",
    "ΔPolar vdW",
    "ΔFavorable"
]


# =========================
# Plot
# =========================
fig = plt.figure(figsize=(13, 11))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.28)

# A
axA = fig.add_subplot(gs[0, 0])
imA = draw_heatmap(
    axA, inter_mat, interaction_cols, master["System"].tolist(),
    "Interaction-type counts across six systems", "YlGnBu"
)
panel_label(axA, "A")
cbarA = fig.colorbar(imA, ax=axA, fraction=0.046, pad=0.03)
cbarA.ax.set_ylabel("Count", rotation=270, labelpad=12)

# B
axB = fig.add_subplot(gs[0, 1])
imB = draw_heatmap(
    axB, pep_mat, pep_xlabels, pep_wide.index.tolist(),
    "Peptide residue contribution", "Blues"
)
panel_label(axB, "B")
cbarB = fig.colorbar(imB, ax=axB, fraction=0.046, pad=0.03)
cbarB.ax.set_ylabel("Contact pairs", rotation=270, labelpad=12)

# C
axC = fig.add_subplot(gs[1, 0])
imC = draw_heatmap(
    axC, hot_mat, hot_wide.columns.tolist(), hot_wide.index.tolist(),
    "HLA hotspot residue contacts", "Purples"
)
panel_label(axC, "C")
cbarC = fig.colorbar(imC, ax=axC, fraction=0.046, pad=0.03)
cbarC.ax.set_ylabel("Contact pairs", rotation=270, labelpad=12)

# D
axD = fig.add_subplot(gs[1, 1])
imD = draw_heatmap(
    axD, delta_mat, delta_xlabels, delta_ylabels,
    "WT/APL delta (APL − WT)", "OrRd"
)
panel_label(axD, "D")
cbarD = fig.colorbar(imD, ax=axD, fraction=0.046, pad=0.03)
cbarD.ax.set_ylabel("Delta", rotation=270, labelpad=12)

fig.tight_layout()

png_out = os.path.join(OUT_DIR, "Supplementary_CoCoMaps2_all.png")
pdf_out = os.path.join(OUT_DIR, "Supplementary_CoCoMaps2_all.pdf")
svg_out = os.path.join(OUT_DIR, "Supplementary_CoCoMaps2_all.svg")

fig.savefig(png_out, dpi=300, bbox_inches="tight")
fig.savefig(pdf_out, bbox_inches="tight")
fig.savefig(svg_out, bbox_inches="tight")

plt.close(fig)

print("✅ Supplementary CoCoMaps2 figures generated:")
print(png_out)
print(pdf_out)
print(svg_out)