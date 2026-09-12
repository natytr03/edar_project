import sys
import pandas as pd
import numpy as np
from scipy.stats import percentileofscore

OUTPUT_FILE = "fst_overview_percentile.txt"
sys.stdout = open(OUTPUT_FILE, "w", encoding="utf-8")

VARIANTSET_FILE = r"skripts_outputs\02_population\variantset\EDAR_Variant_Set_complete_final.csv"
FST_FILE = r"data\raw\edar_fst.weir.fst"
fst_df = pd.read_csv(FST_FILE, sep="\t")

fst_df["WEIR_AND_COCKERHAM_FST"] = pd.to_numeric(fst_df["WEIR_AND_COCKERHAM_FST"], errors="coerce")
fst_ref = fst_df["WEIR_AND_COCKERHAM_FST"].dropna()

print(fst_ref.describe())
print("\n")
print("Number of valid FST values:", len(fst_ref))
print("Minimum:", fst_ref.min())
print("Maximum:", fst_ref.max())
print(fst_ref.quantile([0.90, 0.95, 0.99, 0.995]))
print("\n")

fst_lookup = fst_df[["CHROM", "POS", "WEIR_AND_COCKERHAM_FST"]].copy()

fst_lookup = fst_lookup.rename(columns={"CHROM": "chr", "POS": "pos",
                                        "WEIR_AND_COCKERHAM_FST": "FST"})

# filter variant with fst estimate
df = pd.read_csv(VARIANTSET_FILE)
df["FST_percentile"] = df["FST"].apply(lambda x: percentileofscore(fst_ref, x, kind="weak")
if pd.notna(x) else np.nan)

# calculate the percentile
panel_b = df[df["FST"].notna()].copy()
panel_b["FST_percentile"] = panel_b["FST"].apply(lambda x: percentileofscore(
    fst_ref, x, kind="weak"))

#output
print(panel_b[["rsID", "Protein", "Group", "FST",
               "FST_percentile", "phyloP"]].to_string(index=False))

sys.stdout.close()
