import pandas as pd


variants = pd.read_csv("EDAR_Variant_Set_with_FST_new.csv")

#load phyloP bedgraph
phylop = pd.read_csv(
    "edar_phyloP100way.bedgraph",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "phyloP"]
)

phylop["POS_GRCh37"] = phylop["start"] + 1

variants["POS_GRCh37"] = variants["POS_GRCh37"].astype("Int64")
merged = variants.merge(phylop[["POS_GRCh37", "phyloP"]], on="POS_GRCh37", how="left")

merged.to_csv("EDAR_Variant_Set_with_FST_phyloP_new.csv", index=False)
print("Gespeichert: EDAR_Variant_Set_with_FST_phyloP_new.csv")