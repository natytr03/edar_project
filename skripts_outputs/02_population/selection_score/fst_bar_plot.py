import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

df = pd.read_csv(r"skripts_outputs\02_population\variantset\EDAR_Variant_Set_complete_final.csv")
df = df.drop_duplicates(subset="rsID", keep="first").reset_index(drop=True)

df["fst_computed"] = df["FST"].notna()
df["FST_plot"] = pd.to_numeric(df["FST"], errors="coerce").fillna(0.0)
df["label"] = df["Protein"].fillna(df["rsID"])

# rank fst estimates
df = df.sort_values("FST_plot").reset_index(drop=True)

GROUP_COLOR = {"A": "#BBBBBB",
               "B": "#2980B9",
               "C": "#27AE60", }

V370A_COLOR = "#E74C3C"
V370A_EDGE = "#922B21"

fig, ax = plt.subplots(figsize=(14, 5.5))
fig.patch.set_facecolor("#FAFAFA")
ax.set_facecolor("#F8F9FA")

# ranking the bars
for i, row in df.iterrows():
    is_v370a = row["rsID"] == "rs3827760"

    if is_v370a:
        color = V370A_COLOR
        edgecolor = V370A_EDGE
        lw = 2.0
    else:
        color = GROUP_COLOR.get(row["Group"], "#BBBBBB")
        edgecolor = "white"
        lw = 0.5

    # Schraffur for absence in 1000 gp
    hatch = "/////" if not row["fst_computed"] else ""

    ax.bar(i, row["FST_plot"],
           color=color, hatch=hatch,
           edgecolor=edgecolor, linewidth=lw,
           zorder=2, alpha=0.88)

# highlight/annotate V370A
v370a_idx = df[df["rsID"] == "rs3827760"].index[0]
v370a_fst = df.loc[v370a_idx, "FST_plot"]

ax.annotate(
    f"V370A\nFST = {v370a_fst:.3f}\n(EAS vs. NON-EAS)",
    xy=(v370a_idx, v370a_fst),
    xytext=(v370a_idx - 9, 0.85),
    fontsize=8.5, fontweight="bold", color=V370A_COLOR,
    arrowprops=dict(arrowstyle="->", color=V370A_COLOR, lw=1.5),
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
              edgecolor=V370A_COLOR, alpha=0.92, linewidth=1.5)
)

# p.Met107Val annotation (group b outlier)
met_idx = df[df["label"].str.contains("Met107|M107", na=False)].index
if len(met_idx):
    idx = met_idx[0]
    ax.annotate(
        f"Met107Val\nFST = {df.loc[idx, 'FST_plot']:.3f}",
        xy=(idx, df.loc[idx, "FST_plot"]),
        xytext=(idx - 6, 0.28),
        fontsize=7.5, color="#2980B9",
        arrowprops=dict(arrowstyle="->", color="#2980B9", lw=1.2),
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                  edgecolor="#2980B9", alpha=0.85, linewidth=1.2)
    )

# threshold
ax.axhline(0.177, color="#888888", linestyle="--", linewidth=0.9, alpha=0.7, zorder=1)
ax.text(0.5, 0.184, "FST = 0.177", fontsize=7.5, color="#888888", va="bottom")

ax.set_xticks(range(len(df)))
ax.set_xticklabels(df["label"], rotation=90, fontsize=7)
ax.set_ylabel("FST  (EAS vs. NON-EAS)\nPopulation Differentiation", fontsize=10, labelpad=6)
ax.set_xlabel("EDAR Variant", fontsize=10, labelpad=6)
ax.set_ylim(0, 1.0)
ax.set_xlim(-0.8, len(df) - 0.2)

ax.set_title(
    "FST Population Differentiation with Curated EDAR Variant Set\n"
    "Group A (ClinVar pathogenic) · Group B (neutral common) · "
    "Group C (V370A + LD partners)",
    fontsize=11, fontweight="bold", color="#1A252F",
    pad=10, loc="left"
)

ax.spines[["top", "right"]].set_visible(False)
ax.spines[["left", "bottom"]].set_color("#CCCCCC")
ax.tick_params(colors="#555555", length=3)
ax.grid(axis="y", color="#DDDDDD", linewidth=0.6, alpha=0.8, zorder=0)

# differantiate in label again, lighter label no FST estimate
for i, tick in enumerate(ax.get_xticklabels()):
    tick.set_color(
        "#333333" if df.iloc[i]["fst_computed"] else "#AAAAAA"
    )

legend_content = [mpatches.Patch(facecolor="#BBBBBB", edgecolor="#999",
                                 label="Group A: pathogenic (FST = 0, present in 1000G)"),
                  mpatches.Patch(facecolor="#AAAAAA", edgecolor="#999", hatch="///",
                                 label="Group A: absent from 1000G (FST not computable)"),
                  mpatches.Patch(facecolor="#2980B9",
                                 label="Group B: neutral common variants"),
                  mpatches.Patch(facecolor="#27AE60",
                                 label="Group C: LD partners V370A"),
                  mpatches.Patch(facecolor="#E74C3C",
                                 label="V370A: positive selection signal"), ]

ax.legend(handles=legend_content, fontsize=8, loc="upper left",
          framealpha=0.92, fancybox=True, edgecolor="#CCCCCC")

plt.tight_layout()
plt.savefig("FST_ranked_barplot_v2.png", dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
print("Saved.")
