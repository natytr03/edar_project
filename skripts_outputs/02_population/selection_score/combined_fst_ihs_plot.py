import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec

VARIANTSET_FILE = r"variantset\EDAR_Variant_Set_complete_final.csv"
FST_FILE = r"data\raw\edar_fst.weir.fst"
IHS_FILE = r"data\results\edar_eas_5mb_ihs.ihs.out.100bins.norm"

# load fst data
fst = pd.read_csv(FST_FILE, sep="\t", names=["CHROM", "pos", "FST"], skiprows=1)
fst["pos_mb"] = fst["pos"] / 1e6
fst["FST"] = pd.to_numeric(fst["FST"], errors="coerce").clip(lower=0)
fst = fst.dropna(subset=["FST"])

# load ihs data
ihs = pd.read_csv(IHS_FILE, sep="\t", header=None, names=["id", "pos", "freq", "ihh1", "ihh0",
                                                          "ihs_raw", "ihs_norm", "flag"])
ihs["pos_mb"] = ihs["pos"] / 1e6
ihs["ihs_abs"] = ihs["ihs_norm"].abs()

# extract position to place it in the regionwindow
# define the window
EDAR_START = 109450000 / 1e6
EDAR_END = 109540000 / 1e6

var_df = pd.read_csv(VARIANTSET_FILE)
var_df["POS_GRCh37"] = pd.to_numeric(var_df["POS_GRCh37"], errors="coerce")

# get variants needed for plot
# get V370A
v370a_row = var_df[var_df["rsID"] == "rs3827760"].iloc[0]
V370A_POS = int(v370a_row["POS_GRCh37"])

# LD-Partner
ld_df = var_df[(var_df["r2_V370A"].notna()) & (var_df["rsID"] != "rs3827760") &
               (var_df["r2_V370A"] >= 0.85)].copy()
LD_PARTNERS = {}

for _, row in ld_df.iterrows():
    pos = int(row["POS_GRCh37"])
    rsid = row["rsID"]
    r2 = float(row["r2_V370A"])

    LD_PARTNERS[pos] = (
        rsid,
        f"r²={r2:.2f}"
    )

# Group B
group_b_df = var_df[var_df["Group"] == "B"].copy()

GROUP_B = dict(zip(
    group_b_df["POS_GRCh37"].astype(int),
    group_b_df["protein_variant"]
))

assert group_b_df["POS_GRCh37"].notna().all(), \
    "Check group B, POS_GRCh37 missing."
group_b_df["pos_grch37"] = group_b_df["POS_GRCh37"].astype(int)

# V370A data
v370a_fst = fst[fst["pos"] == V370A_POS]["FST"].values
v370a_ihs = ihs[ihs["pos"] == V370A_POS]["ihs_abs"].values
v370a_fst_val = v370a_fst[0] if len(v370a_fst) else 0.803
v370a_ihs_val = v370a_ihs[0] if len(v370a_ihs) else 1.04


# colorgrouping
def fst_color(v):
    if v > 0.5:
        return "#E74C3C"
    elif v > 0.177:
        return "#E67E22"
    return "#AAAAAA"


def ihs_color(v):
    if v > 3.0:
        return "#E74C3C"
    elif v > 2.0:
        return "#E67E22"
    return "#AAAAAA"


fst["color"] = fst["FST"].apply(fst_color)
ihs["color"] = ihs["ihs_abs"].apply(ihs_color)

fig = plt.figure(figsize=(14, 10))
fig.patch.set_facecolor("#FAFAFA")
gs = gridspec.GridSpec(2, 1, figure=fig, hspace=0.12, height_ratios=[1, 1])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharex=ax1)
for ax in [ax1, ax2]:
    ax.set_facecolor("#F8F9FA")


# variant plotting for both plots
def add_common(ax, y_max, show_label=True):
    ax.axvspan(EDAR_START, EDAR_END, alpha=0.08, color="#08519C", zorder=0)
    ax.axvline(V370A_POS / 1e6, color="#E74C3C", linewidth=0.8, linestyle=":", alpha=0.5, zorder=1)

    if show_label:
        ax.text((EDAR_START + EDAR_END) / 2, y_max * 0.96,
                "EDAR", fontsize=9, color="#08519C",
                ha="center", va="top", fontweight="bold")


def add_variants(ax, df_data, val_col, label_offset_y=0.4):
    for pos_int, (rsid, r2_label) in LD_PARTNERS.items():
        row = df_data[df_data["pos"] == pos_int]
        if len(row) == 0:
            continue
        val = row[val_col].values[0]
        ax.scatter(pos_int / 1e6, val, c="#27AE60", s=80, marker="D", edgecolors="#1B5E20", linewidths=1.2, zorder=5)
        ax.annotate(f"{rsid}\n({r2_label})",
                    xy=(pos_int / 1e6, val),
                    xytext=(pos_int / 1e6 - 0.3, val + label_offset_y),
                    fontsize=7.5, color="#1B5E20", ha="center",
                    arrowprops=dict(arrowstyle="->", color="#27AE60", lw=1.0),
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              edgecolor="#27AE60", alpha=0.85), zorder=7)

    for pos_int, label in GROUP_B.items():
        row = df_data[df_data["pos"] == pos_int]
        if len(row) == 0:
            continue
        val = row[val_col].values[0]
        ax.scatter(pos_int / 1e6, val, c="#2980B9", s=65, marker="s", edgecolors="white", linewidths=0.8, zorder=5,
                   alpha=0.9)
        if label == "M107V":
            ax.annotate(f"{label}\n(Group B)",
                        xy=(pos_int / 1e6, val),
                        xytext=(pos_int / 1e6 + 0.3, val + label_offset_y),
                        fontsize=7.5, color="#1A5276", ha="center",
                        arrowprops=dict(arrowstyle="->", color="#2980B9", lw=1.0),
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                                  edgecolor="#2980B9", alpha=0.85), zorder=7)


# FST PLOT
y_max_fst = 1.0
ax1.set_ylim(-0.02, y_max_fst)
add_common(ax1, y_max_fst)

for color, c, s, alpha, z in [("#AAAAAA", "#AAAAAA", 8, 0.45, 2),
                              ("#E67E22", "#E67E22", 20, 0.8, 3),
                              ("#E74C3C", "#E74C3C", 30, 0.9, 4), ]:
    sub = fst[fst["color"] == color]
    ax1.scatter(sub["pos_mb"], sub["FST"], c=c, s=s,
                alpha=alpha, zorder=z, linewidths=0)

# highlight and annotate V370
ax1.scatter(V370A_POS / 1e6, v370a_fst_val, c="#E74C3C", s=350, marker="*",
            edgecolors="#922B21", linewidths=1.5, zorder=6)
ax1.annotate(f"V370A\nFST = {v370a_fst_val:.3f}",
             xy=(V370A_POS / 1e6, v370a_fst_val),
             xytext=(V370A_POS / 1e6 - 0.7, v370a_fst_val - 0.15),
             fontsize=8.5, fontweight="bold", color="#C0392B",
             arrowprops=dict(arrowstyle="->", color="#C0392B", lw=1.3,
                             connectionstyle="arc3,rad=0.2"),
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                       edgecolor="#C0392B", alpha=0.95, linewidth=1.5), zorder=7)

add_variants(ax1, fst, "FST", label_offset_y=0.08)

# threshoald, like in FST ranking
ax1.axhline(0.177, color="gray", linewidth=0.8, linestyle=":", alpha=0.6, zorder=1)
ax1.text(107.05, 0.19, "FST = 0.177", fontsize=7.5, color="gray")

ax1.set_ylabel("FST  (EAS vs. NON-EAS)\nPopulation Differentiation", fontsize=10, labelpad=6)
ax1.set_title("Combined Selection Scan in EDAR Locus  (Chr2: 107–112 Mb, GRCh37)\n"
              "FST (EAS vs NON-EAS) &' iHS (EAS)  |  1000 Genomes Phase 3",
              fontsize=12, fontweight="bold", color="#1A252F", pad=10, loc="left")
ax1.spines[["top", "right"]].set_visible(False)
ax1.spines[["left", "bottom"]].set_color("#CCCCCC")
ax1.grid(axis="y", color="#DDDDDD", linewidth=0.5, alpha=0.7)
ax1.tick_params(labelbottom=False, colors="#555")

leg1 = [mpatches.Patch(color="#AAAAAA", alpha=0.6, label="FST ≤ 0.2"),
        mpatches.Patch(color="#E67E22", label="FST > 0.2"),
        mpatches.Patch(color="#E74C3C", label="FST > 0.5"),
        mlines.Line2D([0], [0], marker="*", color="w", markersize=11,
                      markerfacecolor="#E74C3C", markeredgecolor="#922B21",
                      label="Group C: V370A"),
        mlines.Line2D([0], [0], marker="D", color="w", markersize=7,
                      markerfacecolor="#27AE60", markeredgecolor="#1B5E20",
                      label="Group C: LD partners (r² > 0.85)"),
        mlines.Line2D([0], [0], marker="s", color="w", markersize=7,
                      markerfacecolor="#2980B9", markeredgecolor="white",
                      label="Group B: Neutral common variants"),
        mpatches.Patch(color="#08519C", alpha=0.15, label="EDAR gene"), ]

ax1.legend(handles=leg1, fontsize=8, loc="upper left", framealpha=0.92, fancybox=True, edgecolor="#CCC", ncol=4)

# IHS PLOT
y_max_ihs = ihs["ihs_abs"].max() + 0.8
ax2.set_ylim(0, y_max_ihs)
add_common(ax2, y_max_ihs)

for color, c, s, alpha, z in [("#AAAAAA", "#AAAAAA", 8, 0.45, 2),
                              ("#E67E22", "#E67E22", 20, 0.8, 3),
                              ("#E74C3C", "#E74C3C", 30, 0.9, 4), ]:
    sub = ihs[ihs["color"] == color]
    ax2.scatter(sub["pos_mb"], sub["ihs_abs"], c=c, s=s, alpha=alpha, zorder=z, linewidths=0)

# V370A
ax2.scatter(V370A_POS / 1e6, v370a_ihs_val, c="#E74C3C", s=350, marker="*",
            edgecolors="#922B21", linewidths=1.5, zorder=6)
ax2.annotate(f"V370A\n|iHS| = {v370a_ihs_val:.2f}",
             xy=(V370A_POS / 1e6, v370a_ihs_val),
             xytext=(V370A_POS / 1e6 - 0.7, v370a_ihs_val + 0.7),
             fontsize=8.5, fontweight="bold", color="#C0392B",
             arrowprops=dict(arrowstyle="->", color="#C0392B", lw=1.3,
                             connectionstyle="arc3,rad=-0.2"),
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                       edgecolor="#C0392B", alpha=0.95, linewidth=1.5), zorder=7)

add_variants(ax2, ihs, "ihs_abs", label_offset_y=0.45)

# threshold
ax2.axhline(2.0, color="#E67E22", linewidth=0.9, linestyle="--", alpha=0.65, zorder=1)
ax2.axhline(3.0, color="#E74C3C", linewidth=0.9, linestyle="--", alpha=0.65, zorder=1)
ax2.text(107.05, 2.05, "|iHS| = 2.0 (top 5%)", fontsize=7.5, color="#E67E22", va="bottom")
ax2.text(107.05, 3.05, "|iHS| = 3.0 (top 1%)", fontsize=7.5, color="#E74C3C", va="bottom")

ax2.set_xlabel("Genomic Position in Mb (Chr2, GRCh37)", fontsize=10, labelpad=6)
ax2.set_ylabel("|iHS| \nIntegrated Haplotype Score", fontsize=10, labelpad=6)
ax2.set_xlim(107.0, 112.0)
ax2.spines[["top", "right"]].set_visible(False)
ax2.spines[["left", "bottom"]].set_color("#CCCCCC")
ax2.grid(axis="y", color="#DDDDDD", linewidth=0.5, alpha=0.7)
ax2.tick_params(colors="#555")

leg2 = [mpatches.Patch(color="#AAAAAA", alpha=0.6, label="|iHS| ≤ 2.0"),
        mpatches.Patch(color="#E67E22", label="|iHS| > 2.0  (top 5%)"),
        mpatches.Patch(color="#E74C3C", label="|iHS| > 3.0  (top 1%)"), ]

ax2.legend(handles=leg2, fontsize=8, loc="upper left", framealpha=0.92, fancybox=True, edgecolor="#CCC")

plt.savefig("Combined_FST_iHS_v5.png", dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
print("Saved.")
