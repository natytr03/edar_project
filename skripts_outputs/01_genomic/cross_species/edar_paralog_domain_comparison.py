import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle

proteins = [
    dict(
        name="EDAR", uniprot="Q9UNE0", length=448,
        signal_peptide=(1, 23),
        ecd=(24, 182),
        tm=(183, 203),
        cytoplasmic=(204, 448),
        death_domain=(358, 431),
        death_domain_struct=(344, 435),
        traf_note=None,
        variant=(370, "V370A (rs3827760)"), ),

    dict(
        name="EDA2R / XEDAR", uniprot="Q9HAV5", length=297,
        signal_peptide=None,
        ecd=(1, 136),
        tm=(137, 157),
        cytoplasmic=(158, 297),
        death_domain=None,
        death_domain_struct=None,
        traf_note=[("TRAF6", 252, 260)],
        variant=None, ),

    dict(
        name="TNFRSF19 / TROY", uniprot="Q9NS68", length=423,
        signal_peptide=(1, 23),
        ecd=(24, 170),
        tm=(171, 193),
        cytoplasmic=(194, 423),
        death_domain=None,
        death_domain_struct=None,
        traf_note="TRAF-binding (cytoplasmic tail)",
        variant=None, ), ]

COLOR_SIGNAL = "#B0B0B0"
COLOR_ECD = "#4C72B0"
COLOR_TM = "#333333"
COLOR_CYTO = "#DCDCDC"
COLOR_DD = "#E8873A"
COLOR_DD_EDGE = "#8C4A00"
COLOR_TRAF = "#7A5195"
COLOR_VARIANT = "#C0392B"

MAX_LEN = max(p["length"] for p in proteins)

fig, ax = plt.subplots(figsize=(12, 6.2), dpi=300)

track_h = 0.5
track_gap = 1.6
y_positions = [len(proteins) * track_gap - i * track_gap for i in range(len(proteins))]

for p, y in zip(proteins, y_positions):
    length = p["length"]

    ax.add_patch(Rectangle((0, y - track_h / 2), length, track_h, linewidth=1.0, edgecolor="#888888",
                           facecolor="#F5F5F5", zorder=1))

    if p["signal_peptide"]:
        s, e = p["signal_peptide"]
        ax.add_patch(Rectangle((s, y - track_h / 2), e - s + 1, track_h, linewidth=0.8, edgecolor="#555555",
                               facecolor=COLOR_SIGNAL, zorder=2))

    s, e = p["ecd"]
    ax.add_patch(Rectangle((s, y - track_h / 2), e - s + 1, track_h,
                           linewidth=0.8, edgecolor="#555555",
                           facecolor=COLOR_ECD, zorder=2))

    s, e = p["tm"]
    ax.add_patch(Rectangle((s, y - track_h / 2), e - s + 1, track_h, linewidth=0.8, edgecolor="#555555",
                           facecolor=COLOR_TM, zorder=2))

    s, e = p["cytoplasmic"]
    ax.add_patch(Rectangle((s, y - track_h / 2), e - s + 1, track_h, linewidth=0.8, edgecolor="#555555",
                           facecolor=COLOR_CYTO, zorder=1.5))

    # death domain (EDAR only)
    if p["death_domain"]:
        s, e = p["death_domain"]
        ax.add_patch(Rectangle((s, y - track_h / 2), e - s + 1, track_h, linewidth=1.3, edgecolor=COLOR_DD_EDGE,
                               facecolor=COLOR_DD, zorder=5))
        ss, se = p["death_domain_struct"]
        ax.add_patch(
            Rectangle((ss, y - track_h / 2 - 0.06), se - ss + 1, track_h + 0.12, linewidth=1.0, edgecolor=COLOR_DD_EDGE,
                      facecolor="none", linestyle=(0, (4, 2)), zorder=6))

    # TRAF-binding (XEDAR / TROY only) as reference
    if p["traf_note"]:
        if isinstance(p["traf_note"], list):
            for tag, s, e in p["traf_note"]:
                ax.add_patch(Rectangle((s, y - track_h / 2), e - s + 1, track_h, linewidth=1.0, edgecolor=COLOR_TRAF,
                                       facecolor=COLOR_TRAF, alpha=0.55, zorder=4))
                ax.text(e + 8, y, f"{tag}-binding\n(no DD)", ha="left", va="center",
                        fontsize=7.5, color=COLOR_TRAF, style="italic")
        else:
            ax.text(length + 10, y, p["traf_note"] + "\n(no DD)", ha="left", va="center",
                    fontsize=7.5, color=COLOR_TRAF, style="italic")

    # hightlight v370a variant
    if p["variant"]:
        pos, label = p["variant"]
        ax.plot([pos, pos], [y - track_h / 2 - 0.10, y + track_h / 2 + 0.22], color=COLOR_VARIANT, lw=1.6, zorder=7)
        ax.plot(pos, y + track_h / 2 + 0.22, marker="*", markersize=12, color=COLOR_VARIANT,
                zorder=8, markeredgecolor="white", markeredgewidth=0.4)
        ax.text(pos, y + track_h / 2 + 0.32, label, ha="center", va="bottom", fontsize=8.5, fontweight="bold",
                color=COLOR_VARIANT)

    ax.text(-14, y, f"{p['name']}", ha="right", va="center", fontsize=10, fontweight="bold")
    ax.text(-14, y - 0.28, f"{p['uniprot']}  ({length} aa)", ha="right", va="center", fontsize=7.5, color="#555555")

ax.set_xlim(-150, MAX_LEN + 130)
ax.set_ylim(0.2, len(proteins) * track_gap + 1.0)
ax.set_xticks(range(0, MAX_LEN + 1, 50))
ax.set_xlabel("Amino acid position", fontsize=10)
ax.tick_params(axis="x", labelsize=8)
ax.set_yticks([])
for spine in ["top", "right", "left"]:
    ax.spines[spine].set_visible(False)
ax.spines["bottom"].set_color("#AAAAAA")

ax.set_title("EDAR vs. Paralogs: Domain architecture comparison\n", fontsize=12.5, fontweight="bold", pad=16)

legend_elements = [
    patches.Patch(facecolor=COLOR_SIGNAL, edgecolor="#555555", label="Signal peptide"),
    patches.Patch(facecolor=COLOR_ECD, edgecolor="#555555", label="Extracellular domain (TNFR-Cys repeats)"),
    patches.Patch(facecolor=COLOR_TM, edgecolor="#555555", label="Transmembrane domain"),
    patches.Patch(facecolor=COLOR_CYTO, edgecolor="#555555", label="Intracellular domain"),
    patches.Patch(facecolor=COLOR_DD, edgecolor=COLOR_DD_EDGE,
                  label="Death Domain (EDAR only; 358–431 Pfam/UniProt/InterPro)"),
    patches.Patch(facecolor="none", edgecolor=COLOR_DD_EDGE, linestyle=(0, (4, 2)),
                  label="Death Domain, structure-based (344–435, CATH/SCOP)"),
    patches.Patch(facecolor=COLOR_TRAF, edgecolor=COLOR_TRAF, alpha=0.55, label="TRAF-binding motif (XEDAR/TROY)"), ]

ax.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.42, -0.14), ncol=2, fontsize=7.8,
          frameon=False)

plt.tight_layout()
plt.savefig("edar_paralog_domain_comparison.png", dpi=300, bbox_inches="tight")
