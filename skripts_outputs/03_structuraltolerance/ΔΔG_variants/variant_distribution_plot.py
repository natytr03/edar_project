import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
import pandas as pd

df = pd.read_csv("EDAR_DD_all_variants_combined.csv")

pathogenic_variants = {"I431T", "R420Q", "C428R", "W434G", "R358Q",
                       "A378T", "G382S", "G405A", "G405D", "L427S"}

common_variants = {"S380R", "T346M"}

df_pos370 = df[df["position"] == 370]
df_pathogenic = df[df["variant"].isin(pathogenic_variants)]
df_common = df[df["variant"].isin(common_variants)]

# FoldX
foldx_pos370 = dict(zip(df_pos370["variant"], df_pos370["foldx_ddG"]))
foldx_pathogenic = dict(zip(df_pathogenic["variant"], df_pathogenic["foldx_ddG"]))
foldx_common = dict(zip(df_common["variant"], df_common["foldx_ddG"]))

# ThermoMPNN
tmpnn_pos370 = dict(zip(df_pos370["variant"], df_pos370["thermompnn_ddG"]))
tmpnn_pathogenic = dict(zip(df_pathogenic["variant"], df_pathogenic["thermompnn_ddG"]))
tmpnn_common = dict(zip(df_common["variant"], df_common["thermompnn_ddG"]))

# zwischenchecks
assert len(foldx_pos370) == 19
assert len(tmpnn_pos370) == 19
assert len(foldx_pathogenic) == 10
assert len(tmpnn_pathogenic) == 10
assert len(foldx_common) == 2
assert len(tmpnn_common) == 2


# beeswarm offset
def beeswarm_offsets(values, spread=0.16):
    n = len(values)
    order = np.argsort(values)
    seq = [0]
    k = 1
    while len(seq) < n:
        seq += [k, -k]
        k += 1
    seq = seq[:n]
    max_abs = max(abs(s) for s in seq) if max(abs(s) for s in seq) > 0 else 1
    offsets = np.zeros(n)
    for rank_pos, orig_idx in enumerate(order):
        offsets[orig_idx] = seq[rank_pos] / max_abs * spread
    return offsets


def plot_panel(ax, pos370, pathogenic, common, tool_name, ylim):
    group_x = {"pos370": 0, "path": 1, "common": 2}

    #Position 370 group
    labels370 = list(pos370.keys())
    vals370 = np.array(list(pos370.values()))
    offs370 = beeswarm_offsets(vals370, spread=0.28)
    is_v370a = np.array([l == "V370A" for l in labels370])

    ax.scatter(
        group_x["pos370"] + offs370[~is_v370a],
        vals370[~is_v370a],
        color="#888888", s=42, alpha=0.8, zorder=3, edgecolors='none')

    ax.scatter(
        group_x["pos370"] + offs370[is_v370a],
        vals370[is_v370a],
        facecolors='#1f77b4', edgecolors='#1f77b4', s=130, zorder=6,
        marker='o', label='V370A')

    ax.scatter(
        group_x["pos370"] + offs370[is_v370a],
        vals370[is_v370a],
        facecolors='none', edgecolors='white', s=130, linewidths=1.2, zorder=7)

    #ClinVar variants group
    labelsP = list(pathogenic.keys())
    valsP = np.array(list(pathogenic.values()))
    offsP = beeswarm_offsets(valsP, spread=0.22)
    xP = group_x["path"] + offsP

    ax.scatter(xP, valsP, color="#d62728", s=55, zorder=3, edgecolors='none')

    clinvar_labels = {"R358Q", "W434G", "I431T", "R420Q", "C428R"}

    clinvar_texts = []
    for x, y, lab in zip(xP, valsP, labelsP):
        if lab in clinvar_labels:
            clinvar_texts.append(
                ax.text(x, y, lab, fontsize=6.5, color="#a01414", va="center",zorder=8))

    # Common variants group
    labelsN = list(common.keys())
    valsN = np.array(list(common.values()))
    offsN = beeswarm_offsets(valsN, spread=0.12)
    xN = group_x["common"] + offsN

    ax.scatter(xN, valsN, color="#2ca02c", s=80, zorder=4, marker='o')

    common_texts = []
    for x, y, lab in zip(xN, valsN, labelsN):
        common_texts.append(
            ax.text(x, y, lab, fontsize=7, color="#1c6e1c", va="center",
                    zorder=8))

    ax.axhline(0, color='black', linewidth=0.7, linestyle='--', alpha=0.6, zorder=1)

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["V370X\n(n=19)", "ClinVar\n(n=10)", "Common variants\n(n=2)"])
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylim(*ylim)
    ax.set_title(tool_name, fontsize=12, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)

    # adjustText, so labels do not overlap
    all_texts = clinvar_texts + common_texts
    if all_texts:
        adjust_text(
            all_texts,
            ax=ax,
            x=list(xP) + list(xN),
            y=list(valsP) + list(valsN),
            arrowprops=dict(arrowstyle='-', color='#888888', lw=0.5, alpha=0.7),
            expand_points=(1.4, 1.6),
            expand_text=(1.2, 1.3),
            force_text=(0.3, 0.5),
            only_move={'points': 'y', 'text': 'xy'})

    return clinvar_texts, common_texts



all_vals = (list(foldx_pos370.values()) + list(foldx_pathogenic.values()) + list(foldx_common.values()) +
            list(tmpnn_pos370.values()) + list(tmpnn_pathogenic.values()) + list(tmpnn_common.values()))
ylim = (-3, 6)

fig, axes = plt.subplots(1, 2, figsize=(10.5, 5), sharey=True)
plot_panel(axes[0], foldx_pos370, foldx_pathogenic, foldx_common, "FoldX", ylim)
plot_panel(axes[1], tmpnn_pos370, tmpnn_pathogenic, tmpnn_common, "ThermoMPNN", ylim)

axes[0].set_ylabel("ΔΔG (kcal/mol)")
axes[0].legend(loc='upper left', fontsize=8, frameon=False)

fig.suptitle("Predicted Stability Effects in the EDAR Death Domain", fontsize=13, y=1.03)

plt.tight_layout()
plt.savefig("variant_distribution.png", dpi=300, bbox_inches='tight')


