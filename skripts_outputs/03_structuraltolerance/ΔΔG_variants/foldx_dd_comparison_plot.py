import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


df = pd.read_csv("EDAR_DD_all_variants_combined.csv")
# only these variants to plot
variants = ["V370A", "R358Q", "W434G", "I431T", "R420Q", "C428R"]

df_selected = df[df["variant"].isin(variants)]
df_selected = df_selected.set_index("variant").loc[variants]
ddg = df_selected["foldx_ddG"].tolist()
sd = df_selected["foldx_sd"].tolist()


# green = tolerable, yellow = discordance, red = confirm pathogenic
colors = ["#2F5597", "#4A7C8C", "#4A7C8C", "#C0392B", "#C0392B", "#C0392B"]
labels_extra = ["positive control\n(V370A)",
                "ambiguous\n(AM 0.42)",
                "discordance\n(ClinVar pathogenic /\n AM benign)",
                "pathogenic\n(confirmed)",
                "pathogenic\n(confirmed)",
                "pathogenic\n(confirmed)",]

fig, ax = plt.subplots(figsize=(9, 5.5))
x = np.arange(len(variants))
ax.bar(x, ddg, yerr=sd, capsize=4, color=colors, edgecolor="black", linewidth=0.6, width=0.6)

ax.axhline(0, color="black", linewidth=0.8)
ax.axhline(2, color="grey", linestyle="--", linewidth=1)
ax.text(len(variants) - 0.4, 2.15, "destabilization\n threshold (~+2 kcal/mol)", fontsize=8, color="grey")

ax.set_ylabel("ΔΔG (kcal/mol)\n(FoldX BuildModel, n=5 repeated runs)", fontsize=11)
ax.set_xticks(x)
ax.set_xticklabels(variants, fontsize=11)
ax.set_title("Folding Stability: V370A vs. Death-Domain-Varianten", fontsize=12, pad=12)

for i, (v, s) in enumerate(zip(ddg, sd)):
    y = v + s + 0.12 if v >= 0 else v - s - 0.12
    va = "bottom" if v >= 0 else "top"
    ax.text(i, y, f"{v:+.2f}", ha="center", va=va,
            fontsize=9, fontweight="bold", color="#222222")

for i, (v, l) in enumerate(zip(ddg, labels_extra)):
    y = v + (sd[i] + 0.55) if v >= 0 else v - (sd[i] + 0.55)
    va = "bottom" if v >= 0 else "top"
    ax.text(i, y, l, ha="center", va=va, fontsize=7.5, color="dimgrey")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(-3, 7.5)

plt.tight_layout()
plt.savefig("foldx_dd_comparison_plot.png", dpi=300)
#print("Saved.")
