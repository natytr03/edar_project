import pandas as pd
import re

variants = pd.read_csv("data/EDAR_Variant_Set_with_FST_phyloP_new.csv")

#load the specific EDAR AlphaMissense table
am = pd.read_csv(
    "data/AlphaMissense_EDAR_Q9UNE0.tsv",
    sep="\t",
    header=None,
    names=["CHROM","POS_GRCh38_AM","REF","ALT","genome",
           "uniprot_id","transcript_id","protein_variant",
           "am_pathogenicity","am_class"]
)

#standardized name notation of variants
aa3to1 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C',
    'Gln':'Q','Glu':'E','Gly':'G','His':'H','Ile':'I',
    'Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P',
    'Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'
}

def convert_protein_notation(p):
    if pd.isna(p):
        return p
    if p.startswith('p.'):
        m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', p)
        if m:
            ref_aa = aa3to1.get(m.group(1), '?')
            pos = m.group(2)
            alt_aa = aa3to1.get(m.group(3), '?')
            return f"{ref_aa}{pos}{alt_aa}"
    return p  # Gruppe A + C bereits im richtigen Format

variants['Protein_AM'] = variants['Protein'].apply(convert_protein_notation)

#validation
changed = variants[variants['Protein'] != variants['Protein_AM']][['rsID','Protein','Protein_AM']]
if len(changed) > 0:
    print("Konvertierte Notation:")
    print(changed.to_string())
    print()

merged = variants.merge(
    am[["protein_variant","am_pathogenicity","am_class"]],
    left_on="Protein_AM",
    right_on="protein_variant",
    how="left"
)

#when there are more, then one AlphaMissense score, choose the highes one
merged = merged.sort_values("am_pathogenicity", ascending=False)
merged = merged.drop_duplicates(subset="rsID", keep="first")
merged = merged.sort_index()

cols = ["rsID","Protein","Consequence","Group","FST","phyloP","am_pathogenicity","am_class"]
merged.to_csv("EDAR_Variant_Set_complete_final.csv", index=False)
print("\nGespeichert: EDAR_Variant_Set_complete.csv")