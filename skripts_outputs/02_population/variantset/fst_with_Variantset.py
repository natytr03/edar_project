import pandas as pd

variants = pd.read_csv("EDAR_Variant_Set_finalset.csv")

fst = pd.read_csv("data/edar_fst.weir.fst", sep="\t")
fst.columns = ["CHR", "POS_GRCh37", "FST"]

variants["POS_GRCh37"] = variants["POS_GRCh37"].astype("Int64")
fst["POS_GRCh37"] = fst["POS_GRCh37"].astype(int)

#merge with fst values
merged = variants.merge(fst[["POS_GRCh37", "FST"]], on="POS_GRCh37", how="left")
merged["FST"] = merged["FST"].clip(lower=0)

merged.to_csv("EDAR_Variant_Set_with_FST_new.csv", index=False)
print("Gespeichert: EDAR_Variant_Set_with_FST_new.csv")