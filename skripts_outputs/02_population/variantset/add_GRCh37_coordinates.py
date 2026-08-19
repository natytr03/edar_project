import pandas as pd
from pyliftover import LiftOver

VARIANTSET_FILE = r"EDAR_Variant_Set_new.csv"
AF_DIR = r"/data/results"
CHAIN_FILE = r"/data/hg38ToHg19.over.chain.gz"
OUT_FINAL_SET = r"EDAR_Variant_Set_finalset.csv"

#change to using liftover instead of offset
def liftover_pos(lo, pos_grch38, chrom="chr2"):
    if pd.isna(pos_grch38):
        return None
    result = lo.convert_coordinate(chrom, int(pos_grch38) - 1)
    if not result:
        return None
    return result[0][1] + 1

#Add GRCh37 coordinates
def add_grch37_coords(variantset, chainfile):
    #OFFSET = 616750  #for EDAR region (chr2)
    lo = LiftOver(chainfile)

    variantset["POS_GRCh38"] = pd.to_numeric(variantset["POS_GRCh38"], errors="coerce")

    #check variant in Group A + B without GRCh37 pos
    grch37_filter = (variantset["Group"].isin(["A", "B"]) &
            variantset["POS_GRCh38"].notna() &
            variantset["POS_GRCh37"].isna())
    #variantset.loc[grch37_filter, "POS_GRCh37"] = variantset.loc[grch37_filter, "POS_GRCh38"] + OFFSET
    variantset.loc[grch37_filter, "POS_GRCh37"] = variantset.loc[grch37_filter, "POS_GRCh38"].apply(lambda pos: liftover_pos(lo, pos))

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

    af_cols = [f"AF_{p}" for p in POPS]
    variantset["af_source"] = pd.NA

    #Group B already set AF
    group_b = variantset["Group"] == "B"
    variantset.loc[group_b, "af_source"] = "gnomAD_selection_criterion"

    eligible = variantset["Group"] != "B"

    #overwrite if column is NaN
    for pop in POPS:
        col = f"AF_{pop}"

        #we choose 1000G value and manual/gnomAD values before 1000G got applied
        old_val = variantset[col].copy()
        mapped_1000g = variantset["POS_GRCh37"].map(af_all[pop])
        apply_1000g = eligible & mapped_1000g.notna()

        #track where the values come from
        variantset.loc[apply_1000g, col] = mapped_1000g[apply_1000g]

        check_col = f"AF_{pop}_1000G_check"
        variantset[check_col] = mapped_1000g

        newly_1000g = apply_1000g & variantset["af_source"].isna()
        variantset.loc[newly_1000g, "af_source"] = "1000G_Phase3_plink2"

        not_eligible = variantset["Group"] == "B"
        no_1000g_value = mapped_1000g.isna()

        no_1000g_but_had_value = (not_eligible | no_1000g_value) & old_val.notna() & variantset["af_source"].isna()
        variantset.loc[no_1000g_but_had_value, "af_source"] = "gnomAD_absent_from_1000G"

    #Group A variants are absent from 1000G, bc they are rare pathogen -> so set AF=0
    variantset[af_cols] = variantset[af_cols].apply(pd.to_numeric, errors="coerce")
    group_a = variantset["Group"] == "A"
    variantset.loc[group_a, af_cols] = variantset.loc[group_a, af_cols].fillna(0)
    variantset.loc[group_a, "af_source"] = "assumed_absent_1000G"

    #get if an af source doesn't have one
    still_missing = variantset["af_source"].isna()
    variantset.loc[still_missing, "af_source"] = "no_AF_available"

    check_cols = [f"AF_{p}_1000G_check" for p in POPS]
    for pop in POPS:
        col, check_col, diff_col = f"AF_{pop}", f"AF_{pop}_1000G_check", f"AF_{pop}_diff"
        variantset.loc[group_b, diff_col] = (
            variantset.loc[group_b, col] - variantset.loc[group_b, check_col]
        ).abs()

    print(f"\naf_source breakdown:\n{variantset['af_source'].value_counts(dropna=False)}")

    diff_cols = [f"AF_{p}_diff" for p in POPS]
    print(f"\nGroup B: max AF-Differenz (gnomAD vs. 1000G) pro Population:")
    print(variantset.loc[group_b, ["rsID"] + diff_cols].to_string(index=False))

    return variantset, af_cols + check_cols + diff_cols

if __name__ == "__main__":

    variantset = pd.read_csv(VARIANTSET_FILE)
    #print(f"Loaded: {len(variantset)} variants")

    variantset = add_grch37_coords(variantset, CHAIN_FILE)
    variantset, af_cols = add_af(variantset, AF_DIR)

    variantset.to_csv(OUT_FINAL_SET, index=False)
    print(f"\nSaved as {OUT_FINAL_SET}")

    lo = LiftOver(CHAIN_FILE)
    result = lo.convert_coordinate("chr2", 108897145 - 1)  # korrekte GRCh38-Position, 0-based
    print(result)