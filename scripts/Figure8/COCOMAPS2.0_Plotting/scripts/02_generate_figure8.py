#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
02_generate_Figure8_COCOMAPS2_maintext_final.py

Generate publication-quality Figure 8 for the revised manuscript.

Figure 8 panels:
A. Interface area
B. Favorable interface interactions
C. Anchor residue interaction contribution
D. HLA groove hotspot remodeling

COCOMAPS 2.0 terminology is used throughout.

Outputs:
figures/Figure8/
    Figure8_COCOMAPS2_maintext_v2.png
    Figure8_COCOMAPS2_maintext_v2.tiff
    Figure8_COCOMAPS2_maintext_v2.pdf
    Figure8_COCOMAPS2_maintext_v2.svg
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Patch


# ============================================================
# 1. Paths
# ============================================================

ROOT = "/ST_PRESICION/dongdian/AIDD-3D-Pred/cocomaps2/"
SUMMARY_DIR = os.path.join(ROOT, "summary")
OUT_DIR = os.path.join(ROOT, "figures", "Figure8")

os.makedirs(OUT_DIR, exist_ok=True)


# ============================================================
# 2. Figure 7-consistent palette
# ============================================================

COLORS = {
    "E6wt":   "#1A4F8A",
    "E6apl":  "#8B3A0F",
    "E7wt1":  "#5BA3D0",
    "E7apl1": "#E08A3C",
    "E7wt2":  "#3A7EBF",
    "E7apl2": "#C1622A"
}

SYSTEM_ORDER = [
    "E6wt",
    "E6apl",
    "E7wt1",
    "E7apl1",
    "E7wt2",
    "E7apl2"
]

PAIR_LABELS = [
    "E6wt / E6apl",
    "E7wt1 / E7apl1",
    "E7wt2 / E7apl2"
]


# ============================================================
# 3. Publication-quality global style
# ============================================================

# Arial first; Linux environments without Arial can fall back safely.
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = [
    "Arial",
    "Liberation Sans",
    "DejaVu Sans"
]

rcParams["font.size"] = 11

rcParams["axes.titlesize"] = 13
rcParams["axes.titleweight"] = "bold"

rcParams["axes.labelsize"] = 12
rcParams["axes.labelweight"] = "bold"

rcParams["xtick.labelsize"] = 10.5
rcParams["ytick.labelsize"] = 10.5

rcParams["axes.linewidth"] = 1.2

rcParams["xtick.major.width"] = 1.1
rcParams["ytick.major.width"] = 1.1
rcParams["xtick.major.size"] = 5
rcParams["ytick.major.size"] = 5

rcParams["legend.fontsize"] = 10

# Keep text editable/vector-friendly in PDF/PS
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42

# SVG text retained as text
rcParams["svg.fonttype"] = "none"


# ============================================================
# 4. Load data
# ============================================================

master_file = os.path.join(
    SUMMARY_DIR,
    "cocomaps_master_summary.csv"
)

anchor_file = os.path.join(
    SUMMARY_DIR,
    "cocomaps_anchor_contact.csv"
)

hotspot_file = os.path.join(
    SUMMARY_DIR,
    "cocomaps_hla_hotspots.csv"
)


for f in [master_file, anchor_file, hotspot_file]:
    if not os.path.exists(f):
        raise FileNotFoundError(
            f"Required input file not found:\n{f}"
        )


master = pd.read_csv(master_file)
anchor = pd.read_csv(anchor_file)
hotspot = pd.read_csv(hotspot_file)


master["System"] = pd.Categorical(
    master["System"],
    categories=SYSTEM_ORDER,
    ordered=True
)

master = (
    master
    .sort_values("System")
    .reset_index(drop=True)
)


anchor["System"] = pd.Categorical(
    anchor["System"],
    categories=SYSTEM_ORDER,
    ordered=True
)

anchor = (
    anchor
    .sort_values(["System", "AnchorPosition"])
    .reset_index(drop=True)
)


hotspot["System"] = pd.Categorical(
    hotspot["System"],
    categories=SYSTEM_ORDER,
    ordered=True
)

hotspot = (
    hotspot
    .sort_values(["System", "HLAResidueNumber"])
    .reset_index(drop=True)
)


# ============================================================
# 5. Helper functions
# ============================================================

def panel_label(ax, label):
    """
    Add publication-style panel letter.
    """
    ax.text(
        -0.14,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=17,
        fontweight="bold",
        va="bottom",
        ha="left",
        clip_on=False
    )


def draw_box_spines(ax):
    """
    Keep all four spines visible and consistent.
    """
    for side in ["top", "right", "bottom", "left"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(1.2)
        ax.spines[side].set_color("black")


def add_value_labels(
    ax,
    bars,
    fmt="{:.1f}",
    fontsize=10,
    offset_fraction=0.012
):
    """
    Add values above bars with scale-aware offset.
    """
    ymin, ymax = ax.get_ylim()
    offset = (ymax - ymin) * offset_fraction

    for b in bars:
        h = b.get_height()

        ax.text(
            b.get_x() + b.get_width() / 2,
            h + offset,
            fmt.format(h),
            ha="center",
            va="bottom",
            fontsize=fontsize,
            fontweight="normal",
            color="#222222"
        )


def heatmap_text_color(
    value,
    vmin,
    vmax
):
    """
    Use white text on dark heatmap cells and dark text
    on light cells for improved readability.
    """
    if vmax <= vmin:
        return "black"

    normalized = (value - vmin) / (vmax - vmin)

    return "white" if normalized >= 0.62 else "#222222"


# ============================================================
# 6. Prepare Panel A/B data
# ============================================================

bar_systems = [
    "E6wt",
    "E6apl",
    "E7wt1",
    "E7apl1",
    "E7wt2",
    "E7apl2"
]

# More compact and symmetric pair layout
bar_positions = [
    0.00, 0.78,
    2.35, 3.13,
    4.70, 5.48
]

pair_centers = [
    0.39,
    2.74,
    5.09
]


interface_values = [
    master.loc[
        master["System"] == x,
        "InterfaceArea_A2"
    ].iloc[0]
    for x in bar_systems
]


favorable_values = [
    master.loc[
        master["System"] == x,
        "FavorableContacts"
    ].iloc[0]
    for x in bar_systems
]


# ============================================================
# 7. Prepare Panel C: anchor heatmap
# ============================================================

anchor_order = [
    "P1",
    "P2",
    "C-term"
]

# Accept original CSV label "Cterm" but display "C-term"
anchor_tmp = anchor.copy()

anchor_tmp["AnchorPosition"] = (
    anchor_tmp["AnchorPosition"]
    .replace({
        "Cterm": "C-term",
        "C-terminal": "C-term",
        "C_terminal": "C-term"
    })
)


anchor_mat = (
    anchor_tmp
    .pivot_table(
        index="System",
        columns="AnchorPosition",
        values="ContactPairs",
        aggfunc="sum",
        observed=False
    )
    .reindex(SYSTEM_ORDER)
    .reindex(columns=anchor_order)
)


anchor_data = anchor_mat.to_numpy(dtype=float)


# ============================================================
# 8. Prepare Panel D: HLA hotspot heatmap
# ============================================================

top_hotspots = (
    hotspot
    .groupby(
        "HLAResidue",
        as_index=False,
        observed=False
    )["ContactPairs"]
    .sum()
    .sort_values(
        "ContactPairs",
        ascending=False
    )
    .head(10)["HLAResidue"]
    .tolist()
)


hotspot_mat = (
    hotspot
    .pivot_table(
        index="System",
        columns="HLAResidue",
        values="ContactPairs",
        aggfunc="sum",
        observed=False
    )
    .reindex(SYSTEM_ORDER)
    .reindex(columns=top_hotspots)
)


hotspot_data = hotspot_mat.to_numpy(dtype=float)


# ============================================================
# 9. Create figure
# ============================================================

fig = plt.figure(
    figsize=(13.2, 9.4)
)

gs = fig.add_gridspec(
    2,
    2,
    height_ratios=[1.00, 1.05],
    width_ratios=[1.00, 1.08],
    hspace=0.34,
    wspace=0.27
)


# ============================================================
# 10. Panel A — Interface area
# ============================================================

axA = fig.add_subplot(gs[0, 0])


barsA = axA.bar(
    bar_positions,
    interface_values,
    color=[COLORS[x] for x in bar_systems],
    edgecolor="black",
    linewidth=0.9,
    width=0.64
)


axA.set_xticks(pair_centers)

axA.set_xticklabels(
    PAIR_LABELS,
    rotation=0,
    fontweight="bold"
)


axA.set_ylabel(
    "Interface area (Å²)"
)

axA.set_title(
    "Interface area",
    pad=9
)


axA.set_ylim(
    0,
    max(interface_values) * 1.13
)


draw_box_spines(axA)
panel_label(axA, "A")


add_value_labels(
    axA,
    barsA,
    fmt="{:.1f}",
    fontsize=10
)


# ============================================================
# 11. Panel B — Favorable interface interactions
# ============================================================

axB = fig.add_subplot(gs[0, 1])


barsB = axB.bar(
    bar_positions,
    favorable_values,
    color=[COLORS[x] for x in bar_systems],
    edgecolor="black",
    linewidth=0.9,
    width=0.64
)


axB.set_xticks(pair_centers)

axB.set_xticklabels(
    PAIR_LABELS,
    rotation=0,
    fontweight="bold"
)


axB.set_ylabel(
    "Favorable contact count"
)

axB.set_title(
    "Favorable interface interactions",
    pad=9
)


axB.set_ylim(
    0,
    max(favorable_values) * 1.14
)


draw_box_spines(axB)
panel_label(axB, "B")


add_value_labels(
    axB,
    barsB,
    fmt="{:.0f}",
    fontsize=10
)


# ============================================================
# 12. Panel C — Anchor residue interaction contribution
# ============================================================

axC = fig.add_subplot(gs[1, 0])


anchor_vmin = np.nanmin(anchor_data)
anchor_vmax = np.nanmax(anchor_data)


imC = axC.imshow(
    anchor_data,
    aspect="auto",
    cmap="Blues",
    vmin=anchor_vmin,
    vmax=anchor_vmax,
    interpolation="nearest"
)


axC.set_xticks(
    np.arange(len(anchor_order))
)

axC.set_xticklabels(
    anchor_order,
    fontweight="bold"
)


axC.set_yticks(
    np.arange(len(SYSTEM_ORDER))
)

axC.set_yticklabels(
    SYSTEM_ORDER
)


axC.set_title(
    "Anchor residue interaction contribution",
    pad=9
)


draw_box_spines(axC)
panel_label(axC, "C")


# Heatmap cell values
for i in range(anchor_data.shape[0]):
    for j in range(anchor_data.shape[1]):

        val = anchor_data[i, j]

        if not np.isnan(val):

            color = heatmap_text_color(
                val,
                anchor_vmin,
                anchor_vmax
            )

            axC.text(
                j,
                i,
                f"{val:.0f}",
                ha="center",
                va="center",
                fontsize=10,
                fontweight="bold",
                color=color
            )


cbarC = fig.colorbar(
    imC,
    ax=axC,
    fraction=0.046,
    pad=0.035
)


cbarC.ax.tick_params(
    labelsize=9.5,
    width=1.0,
    length=4
)

cbarC.outline.set_linewidth(1.0)

cbarC.ax.set_ylabel(
    "Contact pairs",
    rotation=270,
    labelpad=17,
    fontsize=11,
    fontweight="bold"
)


# ============================================================
# 13. Panel D — HLA groove hotspot remodeling
# ============================================================

axD = fig.add_subplot(gs[1, 1])


hotspot_vmin = np.nanmin(hotspot_data)
hotspot_vmax = np.nanmax(hotspot_data)


imD = axD.imshow(
    hotspot_data,
    aspect="auto",
    cmap="Purples",
    vmin=hotspot_vmin,
    vmax=hotspot_vmax,
    interpolation="nearest"
)


axD.set_xticks(
    np.arange(len(top_hotspots))
)

axD.set_xticklabels(
    top_hotspots,
    rotation=45,
    ha="right",
    rotation_mode="anchor"
)


axD.set_yticks(
    np.arange(len(SYSTEM_ORDER))
)

axD.set_yticklabels(
    SYSTEM_ORDER
)


axD.set_title(
    "HLA groove hotspot remodeling",
    pad=9
)


draw_box_spines(axD)
panel_label(axD, "D")


# Heatmap cell values
for i in range(hotspot_data.shape[0]):
    for j in range(hotspot_data.shape[1]):

        val = hotspot_data[i, j]

        if not np.isnan(val):

            color = heatmap_text_color(
                val,
                hotspot_vmin,
                hotspot_vmax
            )

            axD.text(
                j,
                i,
                f"{val:.0f}",
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
                color=color
            )


cbarD = fig.colorbar(
    imD,
    ax=axD,
    fraction=0.046,
    pad=0.035
)


cbarD.ax.tick_params(
    labelsize=9.5,
    width=1.0,
    length=4
)

cbarD.outline.set_linewidth(1.0)

cbarD.ax.set_ylabel(
    "Contact pairs",
    rotation=270,
    labelpad=17,
    fontsize=11,
    fontweight="bold"
)


# ============================================================
# 14. Shared system legend
# ============================================================

legend_handles = [
    Patch(
        facecolor=COLORS[s],
        edgecolor="black",
        linewidth=0.7,
        label=s
    )
    for s in SYSTEM_ORDER
]


fig.legend(
    handles=legend_handles,
    labels=SYSTEM_ORDER,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.008),
    ncol=6,
    frameon=False,
    fontsize=10,
    handlelength=1.3,
    handleheight=0.9,
    columnspacing=1.35
)


# ============================================================
# 15. Layout
# ============================================================

# Reserve bottom area for real legend.
# Removed the previous hexadecimal color-description sentence.
fig.subplots_adjust(
    left=0.075,
    right=0.965,
    top=0.955,
    bottom=0.105,
    hspace=0.36,
    wspace=0.28
)


# ============================================================
# 16. Output paths
# ============================================================

png_out = os.path.join(
    OUT_DIR,
    "Figure8_COCOMAPS2_maintext_final.png"
)

tif_out = os.path.join(
    OUT_DIR,
    "Figure8_COCOMAPS2_maintext_final.tiff"
)

pdf_out = os.path.join(
    OUT_DIR,
    "Figure8_COCOMAPS2_maintext_final.pdf"
)

svg_out = os.path.join(
    OUT_DIR,
    "Figure8_COCOMAPS2_maintext_final.svg"
)


# ============================================================
# 17. Save publication-quality outputs
# ============================================================

# High-resolution raster preview/submission image
fig.savefig(
    png_out,
    dpi=600,
    bbox_inches="tight",
    facecolor="white"
)


# 600-dpi TIFF for journal submission
fig.savefig(
    tif_out,
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
    pil_kwargs={
        "compression": "tiff_lzw"
    }
)


# Vector PDF
fig.savefig(
    pdf_out,
    bbox_inches="tight",
    facecolor="white"
)


# Vector SVG
fig.savefig(
    svg_out,
    bbox_inches="tight",
    facecolor="white"
)


plt.close(fig)


# ============================================================
# 18. Report
# ============================================================

print("=" * 72)
print("✅ Publication-quality Figure 8 generated")
print("   Analysis terminology: COCOMAPS 2.0")
print("=" * 72)

print(f"PNG  (600 dpi): {png_out}")
print(f"TIFF (600 dpi): {tif_out}")
print(f"PDF  (vector):  {pdf_out}")
print(f"SVG  (vector):  {svg_out}")