import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec

# load data
df = pd.read_csv(r"skripts_outputs\02_population\variantset\EDAR_Variant_Set_complete_final.csv")
pop_cols = ["AF_EAS", "AF_AFR", "AF_EUR", "AF_SAS", "AF_AMR"]
df[pop_cols] = df[pop_cols].fillna(0)


# protein for label
def make_label(row):
    protein = str(row.get("Protein", "")).strip()
    rsid = str(row["rsID"]).strip()
    protein = protein.replace("p.", "").replace("nan", "")
    if protein and protein != rsid:
        return f"{protein} · {rsid}"
    return rsid


df["label"] = df.apply(make_label, axis=1)

# sort them after their groups
group_order = {"A": 0, "B": 1, "C": 2}
df["group_rank"] = df["Group"].map(group_order)
df["POS_sort"] = df["POS_GRCh38"].fillna(df["POS_GRCh37"])
df = df.sort_values(["group_rank", "POS_sort"]).reset_index(drop=True)

GROUP_COLORS = {"A": "#C0392B", "B": "#2980B9", "C": "#27AE60"}
POP_LABELS = {"AF_EAS": "EAS", "AF_AFR": "AFR", "AF_EUR": "EUR",
              "AF_SAS": "SAS", "AF_AMR": "AMR"}

cmap_heat = LinearSegmentedColormap.from_list(
    "af", ["#EFF5FB", "#4A90D9", "#1A5276", "#0A1F3A"], N=256)

fig = plt.figure(figsize=(12, 18))
fig.patch.set_facecolor("#FFFFFF")
gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[3.5, 1], hspace=0.4)
ax_heat = fig.add_subplot(gs[0])
ax_bar = fig.add_subplot(gs[1])

# HEATMAP PLOT
mat = df[pop_cols].values
n_var = len(df)
n_pop = len(pop_cols)

im = ax_heat.imshow(mat, aspect="auto", cmap=cmap_heat, vmin=0, vmax=1, interpolation="nearest")

ax_heat.set_xticks(np.arange(n_pop))
ax_heat.set_xticklabels([POP_LABELS[c] for c in pop_cols], fontsize=12, fontweight="bold")
ax_heat.xaxis.set_ticks_position("top")
ax_heat.xaxis.set_label_position("top")
ax_heat.set_xlabel("Population (1000 Genomes Phase 3)", fontsize=11, labelpad=8)
ax_heat.set_yticks(np.arange(n_var))
ax_heat.set_yticklabels(df["label"], fontsize=8)
ax_heat.tick_params(left=False, bottom=False, top=False)

for x in np.arange(-0.5, n_pop, 1):
    ax_heat.axvline(x, color="white", linewidth=1.2)
for y in np.arange(-0.5, n_var, 1):
    ax_heat.axhline(y, color="white", linewidth=0.5)

# highlight the groups
groups = df["Group"].tolist()
prev = None
for i, g in enumerate(groups):
    if g != prev and i > 0:
        ax_heat.axhline(i - 0.5, color="#333333", linewidth=2.0, zorder=5)
    prev = g

# color labels after groups
ax_heat.figure.canvas.draw()
for tick, grp in zip(ax_heat.get_yticklabels(), df["Group"].tolist()):
    tick.set_color(GROUP_COLORS[grp])

for grp, color in GROUP_COLORS.items():
    idxs = [i for i, g in enumerate(groups) if g == grp]
    if not idxs:
        continue
    mid_data = np.mean(idxs)
    top_data = min(idxs) - 0.45
    bot_data = max(idxs) + 0.45
    mid_ax = ax_heat.transData.transform([0, mid_data])
    mid_ax = ax_heat.transAxes.inverted().transform(mid_ax)[1]
    top_ax = ax_heat.transData.transform([0, top_data])
    top_ax = ax_heat.transAxes.inverted().transform(top_ax)[1]
    bot_ax = ax_heat.transData.transform([0, bot_data])
    bot_ax = ax_heat.transAxes.inverted().transform(bot_ax)[1]

    ax_heat.plot([1.02, 1.02], [top_ax, bot_ax], color=color,
                 linewidth=8, solid_capstyle="butt",
                 transform=ax_heat.transAxes, clip_on=False)

    ax_heat.text(1.05, mid_ax, f"Group {grp}",
                 fontsize=9, fontweight="bold", color=color,
                 va="center", ha="left",
                 transform=ax_heat.transAxes, clip_on=False)

# hightlight V370A row
v370a_rows = [i for i, r in enumerate(df["rsID"]) if r == "rs3827760"]
if v370a_rows:
    i = v370a_rows[0]
    for x in range(n_pop):
        ax_heat.add_patch(plt.Rectangle(
            (x - 0.5, i - 0.5), 1, 1,
            fill=False, edgecolor="#E74C3C", linewidth=2.5, zorder=6))
    ax_heat.text(n_pop - 0.5, i + 0.45, "▲ V370A",
                 fontsize=8, color="#E74C3C", fontweight="bold",
                 ha="right", va="top", clip_on=False,
                 transform=ax_heat.transData)

# colorbar for af
cbar = plt.colorbar(im, ax=ax_heat, pad=0.12, shrink=0.6, orientation="vertical", aspect=25)
cbar.set_label("Allele Frequency", fontsize=10)
cbar.ax.tick_params(labelsize=9)

ax_heat.set_title("Allele Frequency Distribution with EDAR Variant Set\n"
                  "(1000 Genomes Phase 3)",
                  fontsize=13, fontweight="bold", pad=30, color="#2C3E50")

# INCLUDING BAR PLOT FOR AF DISTRIBUTION OF V370A
v370a = df[df["rsID"] == "rs3827760"]
if len(v370a):
    row = v370a.iloc[0]
    pops = list(POP_LABELS.values())
    afs = [row[c] for c in pop_cols]
    bar_colors = ["#27AE60" if p in ("EAS", "AMR") else "#AED6F1" for p in pops]

    bars = ax_bar.bar(pops, afs, color=bar_colors, edgecolor="white",
                      linewidth=1.2, width=0.55, zorder=3)
    ax_bar.set_ylim(0, 1.1)
    ax_bar.set_ylabel("Allele Frequency", fontsize=10)
    ax_bar.set_title("V370A (rs3827760): Population-specific Allele Frequency",
                     fontsize=11, fontweight="bold", color="#2C3E50", pad=8)
    ax_bar.axhline(0.05, color="#E74C3C", lw=1.5, ls="--",
                   label="5% threshold", zorder=2)
    ax_bar.set_facecolor("#F8F9FA")
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.grid(axis="y", color="white", lw=1.5, zorder=1)
    ax_bar.tick_params(bottom=False)
    ax_bar.legend(fontsize=9, frameon=False)

    for bar, af in zip(bars, afs):
        ax_bar.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.03,
                    f"{af:.1%}", ha="center", va="bottom",
                    fontsize=10, fontweight="bold",
                    color="#1E8449" if af > 0.5 else "#2C3E50")

legend_content = [mpatches.Patch(color=GROUP_COLORS["A"], label="Group A: ClinVar pathogenic"),
                  mpatches.Patch(color=GROUP_COLORS["B"], label="Group B: Neutral common variants"),
                  mpatches.Patch(color=GROUP_COLORS["C"], label="Group C: V370A &' LD partners"), ]

fig.legend(handles=legend_content, loc="lower center", ncol=3,
           fontsize=9, frameon=True, fancybox=True, framealpha=0.9,
           bbox_to_anchor=(0.5, 0.01))

plt.savefig("AF_Distribution_EDAR_2.png", dpi=180, bbox_inches="tight")
print("Saved.")
