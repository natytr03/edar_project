import sys
import pandas as pd
from scipy.stats import percentileofscore

OUTPUT_FILE = "alphamissense_overview_percentile.txt"
sys.stdout = open(OUTPUT_FILE, "w", encoding="utf-8")

AM_FILE = r"data\raw\AlphaMissense_EDAR_Q9UNE0.tsv"
VARIANTSET_FILE = r"skripts_outputs\02_population\variantset\EDAR_Variant_Set_complete_final.csv"

cols = ["chr", "pos", "ref", "alt", "genome",
        "uniprot", "transcript", "aa_change",
        "am_score", "am_class"]

ref = pd.read_csv(AM_FILE, sep="\t", names=cols)

reference_scores = pd.to_numeric(ref["am_score"], errors="coerce").dropna()

print(ref.shape)
print(ref.head())
print(ref["am_score"].describe())
print("\n")

# from my curated dataset who has alphamissense score to calculate percentile
df = pd.read_csv(VARIANTSET_FILE)
df["am_pathogenicity"] = pd.to_numeric(df["am_pathogenicity"], errors="coerce")
df["phyloP"] = pd.to_numeric(df["phyloP"], errors="coerce")

df_am = df[df["am_pathogenicity"].notna() & df["phyloP"].notna()].copy()

#only for checking
print(df_am[["rsID", "Group", "am_pathogenicity", "am_class", "phyloP"]])
print("n =", len(df_am))
print("\n")

df_am["am_percentile"] = df_am["am_pathogenicity"].apply(
    lambda x: percentileofscore(reference_scores, x, kind="weak"))

# get percentile against reference
reference_scores = pd.to_numeric(ref["am_score"], errors="coerce").dropna()

df_am["am_percentile"] = df_am["am_pathogenicity"].apply(
    lambda x: percentileofscore(reference_scores, x, kind="weak"))

#output
print(df_am[["rsID", "Group", "am_pathogenicity",
             "am_percentile", "am_class", "phyloP"]].to_string(index=False))

sys.stdout.close()
