#!/usr/bin/env python3

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from Bio import SeqIO
from Bio.Seq import Seq

FEL_JSON = r"data\results\species_list_edar_5.final_fel.json"
MEME_JSON = r"data\results\species_list_edar_5.final_meme.json"
CODON_FASTA = r"data\processed\species_list_edar_codon_5.final_aligned.fasta"

HUMAN_KEYWORDS = ("homo_sapiens", "edar", "q9une0")

POS_370_BIO = 370
ZOOM_START_BIO = 360
ZOOM_END_BIO = 380

DOMAINS_PROTEIN = [{"name": "TNFR", "start": 30, "end": 149, "color": "#028090"},
                   {"name": "TM Domain", "start": 185, "end": 219, "color": "#0F3460"},
                   {"name": "Death Domain", "start": 321, "end": 435, "color": "#E63946"}, ]


# get correct position by mapping the protein position to HyPhy alignemnt position
def get_hyphy_site(codon_fasta_path, target_bio_pos, human_keywords=HUMAN_KEYWORDS):
    records = list(SeqIO.parse(codon_fasta_path, "fasta"))
    human = next((r for r in records if any(k in r.id.lower() for k in human_keywords)), None)
    if human is None:
        raise ValueError(f"No Human-Sequence: {[r.id for r in records[:5]]}")

    seq = str(human.seq).upper()
    assert len(seq) % 3 == 0, f"Length {len(seq)} couldn't diveded by 3. No Codon-Alignment!"

    n_codons = len(seq) // 3
    ungapped_protein_pos = 0
    for i in range(n_codons):
        codon = seq[i * 3:i * 3 + 3]
        if codon == "---":
            continue
        if "-" in codon:
            print(f"Check alignemnt at {i + 1}: '{codon}'")
            continue
        ungapped_protein_pos += 1
        if ungapped_protein_pos == target_bio_pos:
            return i + 1

    raise ValueError(f"Position {target_bio_pos} outside sequence ({ungapped_protein_pos})")


# get correct position by mapping the HyPhy alignment position to protein position
def build_hyphy_to_protein_map(codon_fasta_path, human_keywords=HUMAN_KEYWORDS):
    records = list(SeqIO.parse(codon_fasta_path, "fasta"))
    human = next((r for r in records if any(k in r.id.lower() for k in human_keywords)), None)
    if human is None:
        raise ValueError(f"No human sequence: {codon_fasta_path}")

    seq = str(human.seq).upper()
    n_codons = len(seq) // 3
    mapping = {}
    protein_pos = 0
    for i in range(n_codons):
        codon = seq[i * 3:i * 3 + 3]
        site = i + 1
        if codon == "---":
            mapping[site] = (None, "–")
        elif "-" in codon:
            mapping[site] = (None, "?")
        else:
            protein_pos += 1
            aa = str(Seq(codon).translate())
            mapping[site] = (protein_pos, aa)
    return mapping


POS_370_HYPHY = get_hyphy_site(CODON_FASTA, POS_370_BIO)
ZOOM_START = get_hyphy_site(CODON_FASTA, ZOOM_START_BIO)
ZOOM_END = get_hyphy_site(CODON_FASTA, ZOOM_END_BIO)
hyphy_to_protein = build_hyphy_to_protein_map(CODON_FASTA)

# save domains first
DOMAINS = []
for d in DOMAINS_PROTEIN:
    DOMAINS.append({
        "name": d["name"],
        "start": get_hyphy_site(CODON_FASTA, d["start"]),
        "end": get_hyphy_site(CODON_FASTA, d["end"]),
        "color": d["color"],
    })

# FEL DATA FOR PLOTTING
with open(FEL_JSON) as f:
    fel = json.load(f)
fel_content = fel["MLE"]["content"]["0"]
n_fel = len(fel_content)
fel_alpha = np.array([s[0] for s in fel_content])
fel_beta = np.array([s[1] for s in fel_content])
fel_pval = np.array([s[4] for s in fel_content])

with np.errstate(divide='ignore', invalid='ignore'):
    log_ratio = np.where(
        (fel_alpha > 0) & (fel_beta > 0),
        np.log2(fel_beta / fel_alpha),
        np.where(fel_beta == 0, -4, 0)
    )
log_ratio = np.clip(log_ratio, -4, 4)
positions = np.arange(1, n_fel + 1)

# MEME DATA FOR PLOTTING
with open(MEME_JSON) as f:
    meme = json.load(f)
meme_content = meme["MLE"]["content"]["0"]
n_meme = len(meme_content)
meme_pval = np.array([s[6] for s in meme_content])

with np.errstate(divide='ignore', invalid='ignore'):
    meme_neglog = np.where(meme_pval > 0, -np.log10(meme_pval), 10)
meme_neglog = np.clip(meme_neglog, 0, 10)
meme_positions = np.arange(1, n_meme + 1)


def get_fel_colors(log_r, pvals, highlight):
    colors = []
    for i, (lr, pv) in enumerate(zip(log_r, pvals)):
        site = i + 1
        if site == highlight:
            colors.append("#FF6B35")
        elif pv <= 0.05 and lr > 0:
            colors.append("#E63946")
        elif pv <= 0.1 and lr > 0:
            colors.append("#F4A261")
        elif pv <= 0.05 and lr < 0:
            colors.append("#028090")
        elif pv <= 0.1 and lr < 0:
            colors.append("#74C2CE")
        else:
            colors.append("#CBD5E1")
    return colors


def get_meme_colors(pvals, highlight):
    colors = []
    for i, pv in enumerate(pvals):
        site = i + 1
        if site == highlight:
            colors.append("#FF6B35")
        elif pv <= 0.05:
            colors.append("#E63946")
        elif pv <= 0.1:
            colors.append("#F4A261")
        else:
            colors.append("#CBD5E1")
    return colors


colors_fel = get_fel_colors(log_ratio, fel_pval, POS_370_HYPHY)
colors_meme = get_meme_colors(meme_pval, POS_370_HYPHY)


def draw_domain_bar(ax, n_sites_axis, pos370, show_xlabel=True):
    ax.set_facecolor("#F8FAFC")
    ax.set_xlim(0, n_sites_axis + 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    for d in DOMAINS:
        rect = mpatches.FancyBboxPatch(
            (d["start"], 0.15),
            d["end"] - d["start"], 0.7,
            boxstyle="round,pad=0.01",
            facecolor=d["color"], edgecolor="white", linewidth=1.5,
            transform=ax.transData
        )
        ax.add_patch(rect)
        mid = (d["start"] + d["end"]) / 2
        fontsize = 7 if d["name"] == "TM" else 8
        ax.text(mid, 0.5, d["name"],
                ha="center", va="center", fontsize=fontsize,
                color="white", fontweight="bold",
                transform=ax.transData)

    ax.axvline(pos370, color="#7B679A", linewidth=2.0)

    if show_xlabel:
        tick_ax = ax.twiny()
        tick_ax.set_xlim(0, n_sites_axis + 1)
        tick_ax.set_xticks(list(range(50, n_sites_axis + 1, 50)))
        tick_ax.xaxis.set_ticks_position("bottom")
        tick_ax.xaxis.set_label_position("bottom")
        tick_ax.set_xlabel("Human EDAR reference (alignment coordinates)",
                           fontsize=11, color="#1A1A2E", labelpad=14)
        tick_ax.tick_params(axis="x", labelsize=9, pad=2)
        tick_ax.spines[["top", "right", "left"]].set_visible(False)
        tick_ax.spines["bottom"].set_color("#CBD5E1")


# FEL PLOT
fig1 = plt.figure(figsize=(16, 7), facecolor="#F8FAFC")
fig1.suptitle("EDAR Selection Profile: FEL Site-Specific dN/dS", fontsize=14, fontweight="bold", color="#1A1A2E",
              y=1.01)

gs1 = gridspec.GridSpec(2, 1, figure=fig1, height_ratios=[5, 0.4], hspace=0.05)
ax_fel = fig1.add_subplot(gs1[0])
ax_dom1 = fig1.add_subplot(gs1[1])

ax_fel.set_facecolor("#F8FAFC")
ax_fel.axhline(0, color="#64748B", linewidth=0.8, zorder=2)
ax_fel.bar(positions, log_ratio, color=colors_fel, width=1.0, linewidth=0, zorder=3, alpha=0.85)
ax_fel.axvline(POS_370_HYPHY, color="#7B679A", linewidth=1.5, linestyle="--", alpha=0.7, zorder=4)

FEL_PVAL_THRESHOLD = 0.05
N_FEL_LABELS = 3

candidate_sites = [site for site in range(1, n_fel + 1)
                   if fel_pval[site - 1] <= FEL_PVAL_THRESHOLD and site != POS_370_HYPHY]

top_fel_sites = sorted(candidate_sites, key=lambda s: fel_pval[s - 1])[:N_FEL_LABELS]

# Legende
LEGEND_X_MIN = n_fel * 0.55
LEGEND_Y_MIN = 2.0

for site in top_fel_sites:
    y_val = log_ratio[site - 1]
    in_legend_zone = (site >= LEGEND_X_MIN) and (y_val >= LEGEND_Y_MIN)
    offset = -1.0 if in_legend_zone else (0.5 if y_val >= 0 else -0.5)
    ax_fel.annotate(str(site), xy=(site, y_val), xytext=(site + 2, y_val + offset),
                    fontsize=9, color="#1A1A2E", arrowprops=dict(arrowstyle="-", color="#94A3B8", lw=0.8))

ax_fel.annotate(f"Position {POS_370_BIO} (V370A)",
                xy=(POS_370_HYPHY, log_ratio[POS_370_HYPHY - 1]),
                xytext=(POS_370_HYPHY + 20, 2.8),
                fontsize=10, color="#7B679A", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#7B679A", lw=1.2))

ax_fel.set_ylabel("dN/dS (log scale)", fontsize=12, color="#1A1A2E", labelpad=8)
ax_fel.set_xlim(0, n_fel + 1)
ax_fel.set_ylim(-4.8, 4.8)
ax_fel.set_xticks([])
ax_fel.spines[["top", "right"]].set_visible(False)
ax_fel.spines[["left", "bottom"]].set_color("#CBD5E1")

legend_content = [mpatches.Patch(color="#E63946", label="Positive selection (p≤0.05)"),
                  mpatches.Patch(color="#F4A261", label="Positive selection (p≤0.1)"),
                  mpatches.Patch(color="#028090", label="Purifying selection (p≤0.05)"),
                  mpatches.Patch(color="#74C2CE", label="Purifying selection (p≤0.1)"),
                  mpatches.Patch(color="#CBD5E1", label="Neutral"),
                  mpatches.Patch(color="#7B679A", label="Position 370 (V370A)"), ]

ax_fel.legend(handles=legend_content, loc="upper right", fontsize=9, framealpha=0.9, ncol=2)

draw_domain_bar(ax_dom1, n_fel, POS_370_HYPHY, show_xlabel=True)

fig1.savefig("Fig1_FEL_v1.png", dpi=200, bbox_inches="tight", facecolor="#F8FAFC")
print("Saved: Fig1_FEL.png")

# MEME PLOT
fig2 = plt.figure(figsize=(16, 6), facecolor="#F8FAFC")
fig2.suptitle("EDAR Selection Profile: MEME Episodic Positive Selection", fontsize=14, fontweight="bold",
              color="#1A1A2E", y=1.01)

gs2 = gridspec.GridSpec(2, 1, figure=fig2, height_ratios=[5, 0.4], hspace=0.05)
ax_meme = fig2.add_subplot(gs2[0])
ax_dom2 = fig2.add_subplot(gs2[1])

ax_meme.set_facecolor("#F8FAFC")
ax_meme.bar(meme_positions, meme_neglog, color=colors_meme, width=1.0, linewidth=0, alpha=0.85)
ax_meme.axhline(-np.log10(0.1), color="#E63946", linewidth=1.0, linestyle="--", alpha=0.7, label="p = 0.1 threshold")
ax_meme.axvline(POS_370_HYPHY, color="#7B679A", linewidth=1.5, linestyle="--", alpha=0.7)

MEME_PVAL_THRESHOLD = 0.1
N_LABELS = 7

significant_sites = [site for site in range(1, n_meme + 1)
                     if meme_pval[site - 1] <= MEME_PVAL_THRESHOLD and site != POS_370_HYPHY]

top_sites = sorted(significant_sites, key=lambda s: meme_pval[s - 1])[:N_LABELS]

for site in top_sites:
    y_val = meme_neglog[site - 1]
    ax_meme.annotate(str(site), xy=(site, y_val), xytext=(site + 3, y_val + 0.12),
                     fontsize=9, color="#1A1A2E", fontweight="bold")

ax_meme.annotate(f"Position {POS_370_BIO} (V370A)",
                 xy=(POS_370_HYPHY, meme_neglog[POS_370_HYPHY - 1]),
                 xytext=(POS_370_HYPHY + 20, 2.5),
                 fontsize=10, color="#7B679A", fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="#7B679A", lw=1.2))

ax_meme.set_ylabel("dN/dS (log scale)", fontsize=12, color="#1A1A2E", labelpad=8)
ax_meme.set_xlim(0, n_meme + 1)
ax_meme.set_ylim(0, max(meme_neglog) * 1.3 + 0.3)
ax_meme.set_xticks([])
ax_meme.spines[["top", "right"]].set_visible(False)
ax_meme.spines[["left", "bottom"]].set_color("#CBD5E1")
ax_meme.legend(fontsize=9, loc="upper right", framealpha=0.9)

draw_domain_bar(ax_dom2, n_fel, POS_370_HYPHY, show_xlabel=True)

fig2.savefig("Fig2_MEME.png", dpi=200, bbox_inches="tight", facecolor="#F8FAFC")
print("Saved: Fig2_MEME.png")

# SURROUNDING OF POSITION 370
fig3, ax_zoom = plt.subplots(figsize=(10, 6), facecolor="#F8FAFC")
fig3.suptitle(f"EDAR around Position {POS_370_BIO} (V370A)", fontsize=14, fontweight="bold", color="#1A1A2E", y=1.01)

ax_zoom.set_facecolor("#F0F4FF")

zoom_range = range(ZOOM_START - 1, ZOOM_END)
zoom_pos = [i + 1 for i in zoom_range]
zoom_ratio = [log_ratio[i] for i in zoom_range if i < n_fel]
zoom_pval = [fel_pval[i] for i in zoom_range if i < n_fel]
zoom_colors = get_fel_colors(zoom_ratio, zoom_pval, POS_370_HYPHY)

ax_zoom.axhline(0, color="#64748B", linewidth=0.8, zorder=2)
ax_zoom.bar(zoom_pos[:len(zoom_ratio)], zoom_ratio, color=zoom_colors, width=0.7,
            linewidth=0.5, edgecolor="white", zorder=3)
ax_zoom.axvline(POS_370_HYPHY, color="#FF6B35", linewidth=2.0, linestyle="--", alpha=0.8, zorder=4)

for i, pos in enumerate(zoom_pos[:len(zoom_ratio)]):
    protein_pos, aa = hyphy_to_protein.get(pos, (None, None))
    if aa is not None:
        color = "#FF6B35" if pos == POS_370_HYPHY else "#374151"
        weight = "bold" if pos == POS_370_HYPHY else "normal"
        ax_zoom.text(pos, -5.1, aa, ha="center", va="top",
                     fontsize=11, color=color, fontweight=weight, fontfamily="monospace")

ax_zoom.annotate(f"V{POS_370_BIO}A\n(East Asian\nselective sweep)",
                 xy=(POS_370_HYPHY, log_ratio[POS_370_HYPHY - 1]),
                 xytext=(POS_370_HYPHY + 3, 2.8),
                 fontsize=10, color="#FF6B35", fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="#FF6B35", lw=1.5))

ax_zoom.text(ZOOM_START, 3.8, "Purifying selection at this region", fontsize=9, color="#028090", style="italic")

aa_at_sites = zoom_pos[::2]
aa_labels = []
for site in aa_at_sites:
    protein_pos, _ = hyphy_to_protein.get(site, (None, None))
    aa_labels.append(str(protein_pos) if protein_pos is not None else "–")

ax_zoom.set_xticks(aa_at_sites)
ax_zoom.set_xticklabels(aa_labels, fontsize=10, rotation=45)
ax_zoom.set_xlabel("Amino Acid Position", fontsize=12, color="#1A1A2E", labelpad=8)
ax_zoom.set_ylabel("dN/dS (log scale)", fontsize=12, color="#1A1A2E", labelpad=8)
ax_zoom.set_xlim(ZOOM_START - 0.8, ZOOM_END + 0.8)
ax_zoom.set_ylim(-5.5, 4.5)
ax_zoom.spines[["top", "right"]].set_visible(False)
ax_zoom.spines[["left", "bottom"]].set_color("#CBD5E1")

legend_content = [mpatches.Patch(color="#028090", label="Purifying (p≤0.05)"),
                  mpatches.Patch(color="#74C2CE", label="Purifying (p≤0.1)"),
                  mpatches.Patch(color="#CBD5E1", label="Neutral"),
                  mpatches.Patch(color="#FF6B35", label="Position 370 (V370A)"), ]

ax_zoom.legend(handles=legend_content, loc="upper right", fontsize=9, framealpha=0.9)

fig3.tight_layout()
fig3.savefig("Fig3_Zoom_370_v1.png", dpi=200, bbox_inches="tight", facecolor="#F8FAFC")
print("Saved: Fig3_Zoom_370.png")
