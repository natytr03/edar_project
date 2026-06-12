import pandas as pd

VARIANTSET_FILE = r"EDAR_Variant_Set.csv"
AF_DIR = r"data/results"
OUT_FINAL_SET = r"EDAR_Variant_Set_Final.csv"


#Add GRCh37 coordinates
def add_grch37_coords(variantset):
    OFFSET = 616750  #for EDAR region (chr2)

    variantset["POS_GRCh38"] = pd.to_numeric(variantset["POS_GRCh38"], errors="coerce")

    #check variant in Group A + B without GRCh37 pos
    grch37_filter = (variantset["Group"].isin(["A", "B"]) &
            variantset["POS_GRCh38"].notna() &
            variantset["POS_GRCh37"].isna())
    variantset.loc[grch37_filter, "POS_GRCh37"] = variantset.loc[grch37_filter, "POS_GRCh38"] + OFFSET

    #add GRCh37 positions for variants without position or inaccurate
    manual_pos = {
        "rs2470612697": 109513428,
        "rs2470612633": 109513413,
        "rs2470612948": 109513492,
        "rs2470708774": 109539858,
    }
    for rsid, pos in manual_pos.items():
        variantset.loc[variantset["rsID"] == rsid, "POS_GRCh37"] = pos

    variantset.loc[variantset["rsID"] == "rs146567337", "POS_GRCh38"] = 108897116.0
    variantset.loc[variantset["rsID"] == "rs146567337", "POS_GRCh37"] = 109513866.0
    variantset.loc[variantset["rsID"] == "rs146567337", "CHR"] = 2.0

    missing = variantset["POS_GRCh37"].isna().sum()
    print(f"GRCh37 positions missing: {missing}")
    return variantset


#Add AF from 1000G (with PLINK2 output)
def add_af(variantset, af_dir):
    POPS = ["EAS", "AFR", "EUR", "SAS", "AMR"]
    variantset["POS_GRCh37"] = pd.to_numeric(variantset["POS_GRCh37"], errors="coerce")

    af_all = {}
    for pop in POPS:
        af_path = f"{af_dir}/af_{pop}.afreq"
        df = pd.read_csv(af_path, sep="\t")
        df.columns = ["CHROM", "ID", "REF", "ALT", "AF", "OBS_CT"]
        df["POS"] = df["ID"].str.split(":").str[1].astype(float)

        #if a position has more allele choose the most common variant for each position
        df = df.sort_values("AF", ascending=False)
        df = df.drop_duplicates(subset="POS", keep="first")
        af_all[pop] = df.set_index("POS")["AF"]
        #print(f"  {pop}: {len(af_all[pop])} positions loaded")

    #overwrite if column is NaN
    for pop in POPS:
        col = f"AF_{pop}"
        mapped = variantset["POS_GRCh37"].map(af_all[pop])
        variantset[col] = variantset[col].where(variantset[col].notna(), mapped)

    #Group A variants are absent from 1000G, bc they are rare pathogen -> so set AF=0
    af_cols = [f"AF_{p}" for p in POPS]
    variantset.loc[variantset["Group"] == "A", af_cols] = (
        variantset.loc[variantset["Group"] == "A", af_cols].fillna(0)
    )

    return variantset, af_cols

if __name__ == "__main__":

    variantset = pd.read_csv(VARIANTSET_FILE)
    #print(f"Loaded: {len(variantset)} variants")

    variantset = add_grch37_coords(variantset)
    variantset, af_cols = add_af(variantset, AF_DIR)

    variantset.to_csv(OUT_FINAL_SET, index=False)
    print(f"\nSaved as {OUT_FINAL_SET}")