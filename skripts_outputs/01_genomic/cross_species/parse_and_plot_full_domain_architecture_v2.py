import argparse
import re

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle

DD_TARGET = ["Death_EDAR", "Death"]
CRD_TARGET = "TNFR_c6"
CRD_CLUSTER_GAP = 60


def parse_domtblout(path):
    species_data = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(None, 22)
            if len(fields) < 21:
                continue
            target_name = fields[0]
            query_name = fields[3]
            qlen = int(fields[5])
            i_evalue = float(fields[12])
            ali_from, ali_to = int(fields[17]), int(fields[18])
            if query_name not in species_data:
                species_data[query_name] = {"qlen": qlen, "hits": []}
            species_data[query_name]["hits"].append(dict(
                target=target_name, i_evalue=i_evalue,
                ali_from=ali_from, ali_to=ali_to,
            ))
    return species_data


def parse_tm_csv(path):
    tm_data = {}
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split(",", 4)
            if len(parts) < 5:
                continue
            species, start, end, n_helices, note = parts
            tm_data[species] = dict(
                start=int(start) if start else None,
                end=int(end) if end else None,
                n_helices=int(n_helices),
                note=note,
            )
    return tm_data

#in case there are multiple hits
def parse_multi_tm_spans(note):
    m = re.search(r"\(([^)]+)\)", note)
    if not m:
        return []
    spans = []
    for part in m.group(1).split(";"):
        part = part.strip()
        if "-" in part:
            s, e = part.split("-")
            spans.append((int(s), int(e)))
    return spans


def get_death_domain(hits):
    for target in DD_TARGET:
        candidates = [h for h in hits if h["target"] == target]
        if candidates:
            best = min(candidates, key=lambda h: h["i_evalue"])
            return dict(start=best["ali_from"], end=best["ali_to"], i_evalue=best["i_evalue"])
    return None


def get_crd_region(hits):
    candidates = sorted([h for h in hits if h["target"] == CRD_TARGET], key=lambda h: h["ali_from"])
    if not candidates:
        return None
    clusters, current = [], [candidates[0]]
    for h in candidates[1:]:
        if h["ali_from"] - current[-1]["ali_from"] <= CRD_CLUSTER_GAP:
            current.append(h)
        else:
            clusters.append(current)
            current = [h]
    clusters.append(current)
    best_cluster = max(clusters, key=len)
    return dict(start=min(h["ali_from"] for h in best_cluster),
                end=max(h["ali_to"] for h in best_cluster))

#AI-assisted code for formatting names for plot
def prettify_species_name(raw_name):
    m = re.search(r"_e[a]?dar_", raw_name, flags=re.IGNORECASE)
    species_part = raw_name[: m.start()] if m else raw_name
    words = species_part.split("_")
    pretty = " ".join(words)
    return pretty[0].upper() + pretty[1:] if pretty else raw_name


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domtblout", required=True)
    ap.add_argument("--tm-csv", required=True)
    ap.add_argument("--out", default="edar_full_domain_architecture_final")
    ap.add_argument("--species-order", default=None)
    args = ap.parse_args()

    species_data = parse_domtblout(args.domtblout)
    tm_data = parse_tm_csv(args.tm_csv)

    species_list = ([s.strip() for s in args.species_order.split(",")]
                     if args.species_order else list(species_data.keys()))

    parsed = {}
    for sp in species_list:
        hits = species_data[sp]["hits"]
        parsed[sp] = dict(
            qlen=species_data[sp]["qlen"],
            dd=get_death_domain(hits),
            crd=get_crd_region(hits),
            tm=tm_data.get(sp),
            pretty=prettify_species_name(sp),
        )

    print(f"{'species':45s} {'qlen':>5s}  {'CRD':>10s}  {'TM':>12s}  {'DD':>12s}")
    for sp in species_list:
        d = parsed[sp]
        crd_s = f"{d['crd']['start']}-{d['crd']['end']}" if d["crd"] else "none"
        if d["tm"] and d["tm"]["n_helices"] == 1:
            tm_s = f"{d['tm']['start']}-{d['tm']['end']}"
        elif d["tm"]:
            tm_s = f"CHECK({d['tm']['n_helices']}x)"
        else:
            tm_s = "none"
        dd_s = f"{d['dd']['start']}-{d['dd']['end']}" if d["dd"] else "none"
        print(f"{d['pretty']:45s} {d['qlen']:5d}  {crd_s:>10s}  {tm_s:>12s}  {dd_s:>12s}")


    #plot
    n = len(species_list)
    fig_h = max(3.0, 0.55 * n + 1.8)
    fig, ax = plt.subplots(figsize=(12, fig_h), dpi=300)

    track_h = 0.40
    track_gap = 1.0
    y_positions = [n * track_gap - i * track_gap for i in range(n)]
    max_len = max(parsed[s]["qlen"] for s in species_list)

    COLOR_BACKBONE = "#F0F0F0"
    COLOR_DD = "#E8873A"
    COLOR_DD_EDGE = "#8C4A00"
    COLOR_CRD = "#4C72B0"
    COLOR_TM = "#333333"
    COLOR_FLAG = "#C0392B"

    for sp, y in zip(species_list, y_positions):
        d = parsed[sp]
        length = d["qlen"]

        ax.add_patch(Rectangle((0, y - track_h / 2), length, track_h,
                                linewidth=1.0, edgecolor="#888888",
                                facecolor=COLOR_BACKBONE, zorder=1))

        if d["crd"]:
            c = d["crd"]
            ax.add_patch(Rectangle((c["start"], y - track_h / 2), c["end"] - c["start"] + 1, track_h,
                                    linewidth=0.9, edgecolor=COLOR_CRD, facecolor=COLOR_CRD,
                                    alpha=0.45, hatch="////", zorder=2))

        if d["tm"] and d["tm"]["n_helices"] == 1:
            t = d["tm"]
            ax.add_patch(Rectangle((t["start"], y - track_h / 2), t["end"] - t["start"] + 1, track_h,
                                    linewidth=0.9, edgecolor=COLOR_TM, facecolor=COLOR_TM, zorder=3))
        elif d["tm"] and d["tm"]["n_helices"] != 1:
            spans = parse_multi_tm_spans(d["tm"]["note"])
            for s, e in spans:
                ax.add_patch(Rectangle((s, y - track_h / 2), e - s + 1, track_h,
                                        linewidth=1.0, edgecolor=COLOR_FLAG, facecolor="none",
                                        hatch="xxxx", linestyle=(0, (2, 1)), zorder=3))
            ax.text(length + 6, y, f"TM: {d['tm']['n_helices']}x predicted, ambiguous",
                    ha="left", va="center", fontsize=6.5, color=COLOR_FLAG, style="italic")

        if d["dd"]:
            dd = d["dd"]
            ax.add_patch(Rectangle((dd["start"], y - track_h / 2), dd["end"] - dd["start"] + 1, track_h,
                                    linewidth=1.2, edgecolor=COLOR_DD_EDGE, facecolor=COLOR_DD, zorder=4))

        label = d["pretty"]
        ax.text(-14, y, label, ha="right", va="center", fontsize=8,
                fontweight="bold" if "Homo sapiens" in label else "normal", style="italic")

    ax.set_xlim(-190, max_len + 180)
    ax.set_ylim(0.2, n * track_gap + 1.0)
    ax.set_xticks(range(0, int(max_len) + 1, 50))
    ax.set_xlabel("Amino acid position",
                  fontsize=9.5, labelpad=8)
    ax.tick_params(axis="x", labelsize=7.5)
    ax.set_yticks([])
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#AAAAAA")

    ax.set_title("EDAR domain architecture across species\n"
                  "(hmmscan: CRD region + Death Domain and DeepTMHMM: transmembrane helix)",
                  fontsize=12, fontweight="bold", pad=14)

    legend_elements = [
        mpatches.Patch(facecolor=COLOR_CRD, edgecolor=COLOR_CRD, alpha=0.45, hatch="////",
                        label="TNFR-Cys repeat region (PF00020, position-clustered)"),
        mpatches.Patch(facecolor=COLOR_TM, edgecolor=COLOR_TM, label="Transmembrane helix (DeepTMHMM)"),
        mpatches.Patch(facecolor="none", edgecolor=COLOR_FLAG, hatch="xxxx", linestyle=(0, (2, 1)),
                        label="Transmembrane helix, ambiguous prediction (>1 candidate)"),
        mpatches.Patch(facecolor=COLOR_DD, edgecolor=COLOR_DD_EDGE, label="Death Domain (PF24979/PF00531, strict hit)"),
    ]
    ax.legend(handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, -0.14),
              ncol=1, fontsize=8.5, frameon=False)

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.savefig(f"{args.out}_v2.png", dpi=300, bbox_inches="tight")
    print(f"\n[info] saved {args.out}.png")


if __name__ == "__main__":
    main()
