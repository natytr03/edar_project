import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# data from examined species
data = [("Chimpanzee", 0.000, "V", "Primate"),
        ("Gorilla", 0.000, "V", "Primate"),
        ("Elephant", 0.000, "V", "Afrotheria"),
        ("Naked mole-rat", 0.117, "V", "Rodent"),
        ("Alligator", 0.123, "I", "Reptile"),
        ("Ostrich", 0.123, "I", "Bird"),
        ("Bearded dragon", 0.123, "I", "Reptile"),
        ("Python", 0.208, "I", "Reptile"),
        ("Shark", 0.279, "I", "Cartilaginous fish"),
        ("Carp", 0.455, "V", "Bony fish"),
        ("Zebrafish", 0.808, "V", "Bony fish"), ]

names = [d[0] for d in data]
rmsd = [d[1] for d in data]
aa = [d[2] for d in data]
group = [d[3] for d in data]

colors = {"V": "#1565C0", "I": "#EF6C00"}
bar_colors = [colors[a] for a in aa]

fig, ax = plt.subplots(figsize=(8, 6))
y = np.arange(len(names))
ax.barh(y, rmsd, color=bar_colors, edgecolor="black", linewidth=0.6, height=0.65)

for i, (r, g) in enumerate(zip(rmsd, group)):
    ax.text(r + 0.02, i, f"{r:.3f} Å  ({g})",
            va="center", fontsize=8.5, color="dimgrey")

ax.axvline(1.0, color="grey", linestyle="--", linewidth=1)
ax.text(1.05, 0.98, "  \nhighly conserved\n threshold (1 \u00c5)",
        fontsize=8, color="grey", va="top", transform=ax.get_xaxis_transform())

ax.set_yticks(y)
ax.set_yticklabels(names, fontsize=10)
ax.set_xlabel("Backbone RMSD to human (\u00c5)", fontsize=11)
ax.set_title("Structural Conservation of the Death Domain across 11 Species\n"
             "(human = reference, RMSD = 0)", fontsize=12, pad=12)
ax.set_xlim(0, 1.15)
ax.invert_yaxis()

legend_content = [Patch(facecolor=colors["V"], edgecolor="black", label="Valine at position 370"),
                  Patch(facecolor=colors["I"], edgecolor="black", label="Isoleucine at position 370"), ]

ax.legend(handles=legend_content, loc="lower right", bbox_to_anchor=(1, 0.15), fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("cross_species_rmsd.png", dpi=300)
print("Saved.")
