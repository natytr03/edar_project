import argparse
import csv
import math
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.DSSP import DSSP


def parse_domains(domain_args):
    domains = []
    for d in domain_args:
        rng, label = d.split(":")
        start, end = (int(x) for x in rng.split("-"))
        domains.append((start, end, label))

    return domains


def which_domain(pos, domains):
    for start, end, label in domains:
        if start <= pos <= end:
            return f"{label} ({start}-{end})"

    return "outside defined domain"


# check neighboring residues in 5, 8, 10 Å
def neighbor_search_multi_radius(chain, position, radii, seq_near_cutoff=8):
    target_res = chain[position]
    target_atom = target_res["CA"]
    atoms = list(chain.get_atoms())
    ns = NeighborSearch(atoms)

    results = {}
    for radius in radii:
        close_atoms = ns.search(target_atom.coord, radius)
        neighbor_ids = sorted(set(a.get_parent().id[1]
                for a in close_atoms
                if a.get_parent().id[1] != position))

        entries = []
        for nid in neighbor_ids:
            seq_dist = abs(nid - position)

            if seq_dist <= seq_near_cutoff:
                cat = "sequence adjacent"
            else:
                cat = "spatial only"

            entries.append((nid, seq_dist, cat))
        results[radius] = entries

    return results

#Shortest atom distance to residue 370
def shortest_atom_distance(chain, pos_a, pos_b):
    res_a = chain[pos_a]
    res_b = chain[pos_b]
    min_dist = min(atom_a - atom_b
        for atom_a in res_a.get_atoms()
        for atom_b in res_b.get_atoms())

    return min_dist


def ca_distance(chain, pos_a, pos_b):
    ca_a = chain[pos_a]["CA"]
    ca_b = chain[pos_b]["CA"]

    return math.sqrt(sum((x - y) ** 2
                         for x, y in zip(ca_a.coord,ca_b.coord))
                     )

#get the secondary structure and relative solvent
def get_secondstructure_rsa(pdb_file, chain_id, position):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure",pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file, dssp="mkdssp", file_type="PDB")
    key = (chain_id,(" ", position, " "))

    if key in dssp:
        return (dssp[key][2], dssp[key][3])

    return None, None


def format_residue_list(entries):
    if not entries:
        return ""

    return ",".join(str(nid)
        for nid, _, _ in entries)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wt", required=True)
    ap.add_argument("--mut", required=True)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--position", type=int, required=True)
    ap.add_argument("--radii", type=float, nargs="+", default=[5, 8, 10])
    ap.add_argument("--domains", nargs="*", default=[])
    ap.add_argument("--seq-near-cutoff", type=int, default=8)
    ap.add_argument("--out", default="wt_vs_mut_comparison.csv")

    args = ap.parse_args()

    pos = args.position

    #getting all results for csv
    parser = PDBParser(QUIET=True)
    wt_structure = parser.get_structure("wt", args.wt)
    mut_structure = parser.get_structure("mut", args.mut)

    wt_chain = wt_structure[0][args.chain]
    mut_chain = mut_structure[0][args.chain]


    domains = parse_domains(args.domains)
    domain_label = (which_domain(pos, domains)
                    if domains
                    else "(no domain)")

    wt_neighbors = neighbor_search_multi_radius(wt_chain, pos, args.radii, args.seq_near_cutoff)
    mut_neighbors = neighbor_search_multi_radius(mut_chain, pos, args.radii, args.seq_near_cutoff)


    ss_wt, rsa_wt = get_secondstructure_rsa(args.wt, args.chain, pos)
    ss_mut, rsa_mut = get_secondstructure_rsa(args.mut, args.chain, pos)

    #saved all results in csv
    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow(["metric", "WT (Val370)", "V370A (Ala370)", "Difference"])
        writer.writerow(["Domain", domain_label, domain_label,""])

        for radius in args.radii:
            wt_ids = format_residue_list(wt_neighbors[radius])
            mut_ids = format_residue_list(mut_neighbors[radius])

            wt_set = {nid
                      for nid, _, _ in wt_neighbors[radius]
                      }

            mut_set = {nid
                       for nid, _, _ in mut_neighbors[radius]
                       }

            added = sorted(mut_set - wt_set)
            removed = sorted(wt_set - mut_set)

            difference_parts = []
            if added:
                difference_parts.append("added: " + ",".join(str(x)
                                                             for x in added) )

            if removed:
                difference_parts.append("removed: " + ",".join(str(x)
                                                               for x in removed))

            if not difference_parts:
                difference = "None"
            else:
                difference = "; ".join(difference_parts)

            writer.writerow([f"Neighboring residues ({radius:g} Å)", wt_ids, mut_ids,difference])


        if (rsa_wt is not None and rsa_mut is not None):
            rsa_diff = rsa_mut - rsa_wt

            writer.writerow(["Relative solvent accessibility (RSA)", f"{rsa_wt:.2f}", f"{rsa_mut:.2f}", f"{rsa_diff:+.2f}"])
        else:
            writer.writerow(["Relative solvent accessibility (RSA)","","",""])


        ss_wt_out = (ss_wt
                     if ss_wt is not None
                     else "-")
        ss_mut_out = (ss_mut
                      if ss_mut is not None
                      else "-")

        ss_difference = ("None"
                         if ss_wt_out == ss_mut_out
                         else "Changed")

        writer.writerow(["Secondary structure (DSSP)", ss_wt_out, ss_mut_out, ss_difference])


        radius_10 = 10.0
        if radius_10 in wt_neighbors:
            wt_10 = {nid
                     for nid, _, _ in wt_neighbors[radius_10]}
            mut_10 = {nid
                      for nid, _, _ in mut_neighbors[radius_10]}


            radius_8 = 8.0
            if radius_8 in wt_neighbors:
                wt_8 = {nid
                        for nid, _, _ in wt_neighbors[radius_8]}
                mut_8 = {nid
                         for nid, _, _ in mut_neighbors[radius_8]}

                wt_additional = sorted(wt_10 - wt_8)

                mut_additional = sorted(mut_10 - mut_8)

                writer.writerow(["Additional residues examined (10 Å)",",".join(str(x) for x in wt_additional),
                                 ",".join(str(x) for x in mut_additional), "None" if wt_additional == mut_additional
                                 else "Changed"])


        residues_to_compare = sorted(set(nid
                                         for nid, _, _ in wt_neighbors.get(10.0, [])) |
                                     set(nid
                                         for nid, _, _ in mut_neighbors.get(10.0, [])))


        wt_8 = {nid
                for nid, _, _ in wt_neighbors.get(8.0, [])}

        mut_8 = {nid
                 for nid, _, _ in mut_neighbors.get(8.0, [])}

        additional_residues = sorted((set(residues_to_compare) - wt_8 - mut_8))

        additional_residues = sorted((set(nid
                                          for nid, _, _ in wt_neighbors.get(10.0, [])) |
                                      set(nid
                                          for nid, _, _ in mut_neighbors.get(10.0,[]))) -
                                     (wt_8 & mut_8))

        for nid in additional_residues:
            wt_ca = ""
            wt_atom = ""

            mut_ca = ""
            mut_atom = ""

            if nid in wt_chain:
                wt_ca = f"{ca_distance(wt_chain, pos, nid):.2f}"
                wt_atom = f"{shortest_atom_distance(wt_chain, pos, nid):.2f}"

            if nid in mut_chain:
                mut_ca = f"{ca_distance(mut_chain, pos, nid):.2f}"
                mut_atom = f"{shortest_atom_distance(mut_chain, pos, nid):.2f}"

            writer.writerow([f"Distances to residue {nid}", f"Cα={wt_ca}; atom={wt_atom}", f"Cα={mut_ca}; atom={mut_atom}", ""])


if __name__ == "__main__":
    main()
