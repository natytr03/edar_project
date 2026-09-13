import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# parse all values needed for the plot
df = pd.read_csv(r"skripts_outputs\02_population\variantset\EDAR_Variant_Set_complete_final.csv")
df["FST"] = pd.to_numeric(df["FST"], errors="coerce")
df["phyloP"] = pd.to_numeric(df["phyloP"], errors="coerce")
df["am_pathogenicity"] = pd.to_numeric(df["am_pathogenicity"], errors="coerce")

# define for groups
GROUP_STYLE = {"A": dict(color="#C0392B", marker="o", label="Group A: ClinVar pathogenic"),
               "B": dict(color="#2980B9", marker="s", label="Group B: Neutral common"),
               "C": dict(color="#27AE60", marker="o", label="Group C: LD partners"), }

# alphamissense threshold
AM_LIKELY_BENIGN = 0.340
AM_AMBIGUOUS = 0.564

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 13), gridspec_kw={"hspace": 0.45})
fig.patch.set_facecolor("#FAFAFA")

# PLOT PANEL A) phyloP × AlphaMissense
ax1.set_facecolor("#F8F9FA")

# filter only missense variants with am score
miss = df[(df["Consequence"] == "missense_variant") & df["phyloP"].notna() &
          df["am_pathogenicity"].notna()].copy()

# set up for am classifications
ax1.axhspan(0, AM_LIKELY_BENIGN, alpha=0.07, color="#2196F3", zorder=0, label="_nolegend_")
ax1.axhspan(AM_LIKELY_BENIGN, AM_AMBIGUOUS, alpha=0.07, color="#FFC107", zorder=0)
ax1.axhspan(AM_AMBIGUOUS, 1.0, alpha=0.07, color="#F44336", zorder=0)

ax1.text(9.3, 0.17, "likely\nbenign", fontsize=7.5, color="#1565C0", ha="center", style="italic")
ax1.text(9.3, 0.45, "ambiguous", fontsize=7.5, color="#E65100", ha="center", style="italic")
ax1.text(9.3, 0.78, "likely\npathogenic", fontsize=7.5, color="#B71C1C", ha="center", style="italic")

ax1.axhline(AM_LIKELY_BENIGN, color="#888", linewidth=0.7, linestyle="--", alpha=0.5)
ax1.axhline(AM_AMBIGUOUS, color="#888", linewidth=0.7, linestyle="--", alpha=0.5)

# Group A variants
a_grey = df[(df["Group"] == "A") & (df["Consequence"] == "missense_variant") &
            df["phyloP"].notna() & df["am_pathogenicity"].notna()]

ax1.scatter(
    a_grey["phyloP"], a_grey["am_pathogenicity"],
    c="#BBBBBB", s=70, marker="o",
    edgecolors="#999999", linewidths=0.6,
    zorder=2, alpha=0.55,
    label="Group A: missense (no selection signal)"
)

for _, row in a_grey.iterrows():
    ax1.annotate(
        row["Protein"],
        (row["phyloP"], row["am_pathogenicity"]),
        xytext=(3, 4), textcoords="offset points",
        fontsize=6.5, color="#555",
    )

# only for comparison, variants without am as triangle
SEPARATOR_Y = -0.04  # Y-Position der Dreiecke
FLOOR_Y = -0.01  # Y-Position der Trennlinie

a_stop = df[(df["Group"] == "A") &
            (df["Consequence"].isin(["stop_gained", "splice_donor_variant", "splice_acceptor_variant"])) &
            df["phyloP"].notna()]

ax1.scatter(
    a_stop["phyloP"],
    [SEPARATOR_Y] * len(a_stop),
    c="#BBBBBB", s=55, marker="v",
    edgecolors="#999999", linewidths=0.6,
    zorder=2, alpha=0.6, clip_on=False,
    label="Group A: stop/splice (no AM score)"
)

ax1.axhline(FLOOR_Y, color="#AAAAAA", linewidth=1.2, linestyle="-", alpha=0.8, zorder=5)
ax1.text(
    0.01, FLOOR_Y + 0.002,
    "▼  no AM score available (stop / splice variants)",
    transform=ax1.get_yaxis_transform(),
    fontsize=7.5, color="#888888", style="italic", va="bottom"
)

ax1.axhspan(-0.07, FLOOR_Y, color="#EEEEEE", alpha=0.6, zorder=0)
ax1.set_ylim(-0.07, 1.08)

# Group B and C
for grp in ["B", "C"]:
    style = GROUP_STYLE[grp]
    sub = miss[miss["Group"] == grp]
    if len(sub) == 0:
        continue
    ax1.scatter(
        sub["phyloP"], sub["am_pathogenicity"],
        c=style["color"], s=90, marker=style["marker"],
        edgecolors="white", linewidths=0.8,
        zorder=3, alpha=0.88, label=style["label"]
    )
    if grp == "B":
        for _, row in sub.iterrows():
            ax1.annotate(
                row["Protein"],
                (row["phyloP"], row["am_pathogenicity"]),
                xytext=(3, 4), textcoords="offset points",
                fontsize=6.5, color="#555"
            )

# hightlight V370A
v370a = df[df["rsID"] == "rs3827760"]
if len(v370a) and v370a.iloc[0]["am_pathogenicity"] > 0:
    r = v370a.iloc[0]
    ax1.scatter(r["phyloP"], r["am_pathogenicity"],
                c="#E74C3C", s=450, marker="*",
                edgecolors="#922B21", linewidths=1.5, zorder=6)
    ax1.annotate(
        f"V370A\nphyloP = {r['phyloP']:.2f}\nAM = {r['am_pathogenicity']:.3f}",
        xy=(r["phyloP"], r["am_pathogenicity"]),
        xytext=(r["phyloP"] - 2.8, r["am_pathogenicity"] - 0.08),
        fontsize=8.5, fontweight="bold", color="#C0392B",
        arrowprops=dict(arrowstyle="->", color="#C0392B", lw=1.4),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="#C0392B", alpha=0.92, linewidth=1.5)
    )

ax1.set_xlabel("phyloP Conservation Score\n(negative = fast evolving, positive = conserved)", fontsize=10, labelpad=6)
ax1.set_ylabel("AlphaMissense Pathogenicity Score", fontsize=10, labelpad=6)
ax1.set_title("A)  Conservation vs. Predicted Pathogenicity\n"
              "All missense EDAR variants  |  AlphaMissense thresholds shown",
              fontsize=11, fontweight="bold", color="#1A252F", pad=10, loc="left")
ax1.set_ylim(-0.05, 1.08)
ax1.spines[["top", "right"]].set_visible(False)
ax1.grid(color="#DDDDDD", linewidth=0.5, alpha=0.7)
handles, labels = ax1.get_legend_handles_labels()
v370a_handle = mlines.Line2D([0], [0], marker="*", color="w", markersize=13,
                             markerfacecolor="#E74C3C", markeredgecolor="#922B21",
                             label="V370A: positive selection (★)")
handles.append(v370a_handle)
labels.append("V370A: positive selection (★)")
ax1.legend(handles=handles, labels=labels, fontsize=7, framealpha=0.9, fancybox=True,
           edgecolor="#CCC", loc="upper left")

# PLOT PANEL B) FST × phyloP (only variants with estimated FST values)
ax2.set_facecolor("#F8F9FA")
bc = df[df["Group"].isin(["B", "C"]) & df["phyloP"].notna()].copy()

# check in Group A, if there is a variant with an FST estimate
a_hasFST = df[(df["Group"] == "A") & df["FST"].notna() & df["phyloP"].notna()].copy()

ax2.scatter(
    a_hasFST["phyloP"], a_hasFST["FST"],
    c="#CCCCCC", s=60, marker="o",
    edgecolors="#AAAAAA", linewidths=0.6,
    zorder=2, alpha=0.5,
    label="Group A: ClinVar pathogenic"
)

for _, row in a_hasFST.iterrows():
    ax2.annotate(row["Protein"] if pd.notna(row["Protein"]) else "",
                 xy=(row["phyloP"], row["FST"]),
                 xytext=(3, 4), textcoords="offset points",
                 fontsize=6.5, color="#555"
                 )

# Quadranten
ax2.axvspan(-6, 3, alpha=0.04, color="#E67E22", zorder=0)  # variabel
ax2.axvspan(3, 11, alpha=0.04, color="#27AE60", zorder=0)  # konserviert
ax2.axhline(0.15, color="#888", linewidth=0.8, linestyle="--", alpha=0.5)
ax2.axvline(3.0, color="#888", linewidth=0.8, linestyle="--", alpha=0.5)

ax2.text(7.5, 0.88, "conserved\n+ highly differentiated", fontsize=8, color="#1A5E20", ha="center", style="italic")
ax2.text(-1, 0.5, "variable\n+ highly differentiated", fontsize=8, color="#888", ha="center", style="italic")
ax2.text(6, 0.06, "conserved\n+ weakly differentiated", fontsize=8, color="#888", ha="center", style="italic")
ax2.text(-1, 0.06, "variable\n+ weakly differentiated", fontsize=8, color="#888", ha="center", style="italic")

# Group B
b = bc[bc["Group"] == "B"]
for _, row in b.iterrows():
    color = "#2980B9" if pd.isna(row["am_pathogenicity"]) else "#2980B9"
    ax2.scatter(row["phyloP"], row["FST"],
                c=color, s=90, marker="s",
                edgecolors="white", linewidths=0.8,
                zorder=3, alpha=0.85)
    ax2.annotate(row["Protein"] if pd.notna(row["Protein"]) else "",
                 xy=(row["phyloP"], row["FST"]),
                 xytext=(3, 4), textcoords="offset points",
                 fontsize=6.5, color="#555")

# Group C LD-Partner
ld = bc[(bc["Group"] == "C") & (bc["rsID"] != "rs3827760")]
if len(ld):
    ax2.scatter(ld["phyloP"], ld["FST"],
                c="#A8D5A2", s=100, marker="o",
                edgecolors="#27AE60", linewidths=1.2,
                zorder=4, alpha=0.85, label="Group C — LD partners")
    for _, row in ld.iterrows():
        ax2.annotate(
            f"LD (r²={row['r2_V370A']:.2f})" if pd.notna(row.get('r2_V370A')) else "LD partner",
            xy=(row["phyloP"], row["FST"]),
            xytext=(4, 5), textcoords="offset points",
            fontsize=6.5, color="#27AE60", style="italic"
        )

# V370A Stern
if len(v370a):
    r = v370a.iloc[0]
    ax2.scatter(r["phyloP"], r["FST"], c="#E74C3C", s=450, marker="*",
                edgecolors="#922B21", linewidths=1.5, zorder=6)
    ax2.annotate(
        "V370A\nphyloP=8.96\nFST=0.80",
        xy=(r["phyloP"], r["FST"]),
        xytext=(r["phyloP"] - 2.5, r["FST"] - 0.12),
        fontsize=8.5, fontweight="bold", color="#C0392B",
        arrowprops=dict(arrowstyle="->", color="#C0392B", lw=1.4),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="#C0392B", alpha=0.9)
    )

ax2.set_xlabel("phyloP Conservation Score\n(negative = fast evolving, positive = conserved)", fontsize=10, labelpad=6)
ax2.set_ylabel("FST  (EAS vs. NON-EAS)\nPopulation Differentiation", fontsize=10, labelpad=6)
ax2.set_title("B)  Conservation vs. Population Differentiation\n"
              "Group A + Group B + Group C (V370A ★)",
              fontsize=11, fontweight="bold", color="#1A252F", pad=10, loc="left")
ax2.set_ylim(-0.05, 0.95)
ax2.spines[["top", "right"]].set_visible(False)
ax2.grid(color="#DDDDDD", linewidth=0.5, alpha=0.7)

leg2 = [mlines.Line2D([0], [0], marker="o", color="w", markersize=8,
                      markerfacecolor="#CCCCCC", markeredgecolor="#AAAAAA",
                      label="Group A: ClinVar pathogenic (FST measured)"),
        mlines.Line2D([0], [0], marker="s", color="w", markersize=8,
                      markerfacecolor="#2980B9", markeredgecolor="white",
                      label="Group B: Neutral common variants"),
        mlines.Line2D([0], [0], marker="o", color="w", markersize=8,
                      markerfacecolor="#A8D5A2", markeredgecolor="#27AE60",
                      label="Group C: LD partners V370A"),
        mlines.Line2D([0], [0], marker="*", color="w", markersize=10,
                      markerfacecolor="#E74C3C", markeredgecolor="#922B21",
                      label="V370A: positive selection signal"), ]

ax2.legend(handles=leg2, fontsize=6, framealpha=0.9, fancybox=True, edgecolor="#CCC", loc="upper left")

plt.savefig("Two_Correlation_Plots_v7.png", dpi=180, bbox_inches="tight")
print("\nSaved.")
