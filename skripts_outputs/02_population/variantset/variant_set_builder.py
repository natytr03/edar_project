import pandas as pd

CLINVAR_FILE  = r"data\raw\clinvar_variant.txt"
#gnomAD v4.1 (downloaded from gnomAD browser (https://gnomad.broadinstitute.org)
GNOMAD_FILE   = r"data\raw\gnomAD_B_all.csv"
GROUPC_FILE   = r"EDAR_Gruppe_C_new.csv"

OUT_A         = r"EDAR_GroupA_new.csv"
OUT_B         = r"EDAR_GroupB_new.csv"
OUT_VARIANTSET    = r"EDAR_Variant_Set_new.csv"


#GROUP A: Clinically described variants (from ClinVar)
def build_group_a(path):
    df = pd.read_csv(path, sep="\t")
    #print(f"Loaded: {len(df)} variants")

    #filter only variants with rsID
    df = df[df["dbSNP ID"].str.startswith("rs", na=False)]
    #print(f"After rsID filter: {len(df)}")

    #filter only SNVs
    df = df[df["Variant type"].str.contains("single nucleotide", case=False, na=False)]
    #print(f"After SNV filter: {len(df)}")

    #priority in order to choose which variant (nonsense > splice > missense)
    def priority(row):
        cons = str(row["Molecular consequence"]).lower()
        if "nonsense" in cons:   return 1
        elif "splice" in cons:   return 2
        elif "missense" in cons: return 3
        else:                    return 4

    #rank pathogenic
    def path_rank(row):
        c = str(row["Germline classification"]).lower()
        if "pathogenic/likely" in c:   return 1
        elif c == "pathogenic":         return 2
        elif "likely pathogenic" in c:  return 3
        else:                           return 4

    df["priority"]  = df.apply(priority, axis=1)
    df["path_rank"] = df.apply(path_rank, axis=1)
    df["GRCh38Location"] = pd.to_numeric(df["GRCh38Location"], errors="coerce")

    nonsense = df[df["priority"] == 1].sort_values("path_rank")
    splice   = df[df["priority"] == 2].sort_values("path_rank")
    missense = df[df["priority"] == 3].sort_values(["path_rank", "GRCh38Location"])

    df_a = pd.concat([nonsense, splice, missense]).head(30).copy()
    df_a["Group"] = "A"

    #print(f"After selection (head 30): {len(df_a)}")

    #deduplicate by rsID and keep highest pathogenic one
    def path_rank_dedup(c):
        c = str(c).lower()
        if c == "pathogenic": return 1
        elif "pathogenic/likely" in c: return 2
        elif "likely pathogenic" in c: return 3
        else: return 4

    df_a["_rank"] = df_a["Germline classification"].apply(path_rank_dedup)
    df_a = (df_a
        .sort_values("_rank")
        .drop_duplicates(subset="dbSNP ID", keep="first")
        .drop(columns="_rank")
        .reset_index(drop=True))

    #extra removal splice acceptor
    df_a = df_a[df_a["dbSNP ID"] != "rs757233170"]

    #build clean output
    output_cols = {
        "Name":                    "Variant_HGVS",
        "Protein change":          "Protein_Change",
        "dbSNP ID":                "rsID",
        "GRCh38Chromosome":        "CHR",
        "GRCh38Location":          "POS_GRCh38",
        "Molecular consequence":   "Consequence",
        "Germline classification": "ClinVar_Class",
        "Condition(s)":            "Condition",
        "Group":                   "Group",
    }
    df_out = df_a[[c for c in output_cols if c in df_a.columns]].rename(columns=output_cols)
    df_out.to_csv(OUT_A, index=False)

    #print(f"\nGroup A: {len(df_out)} variants saved as {OUT_A}")
    #print(df_out["Consequence"].value_counts().to_string())
    return df_out


#GROUP B: Frequent population variants (with gnomAD AF > 1%)
def build_group_b(path):
    df = pd.read_csv(path)
    #print(f"Loaded: {len(df)} variants")

    #calculate per-population AF (allele frequency) from allele counts
    df["Allele Frequency"] = pd.to_numeric(df["Allele Frequency"], errors="coerce")
    pops = {
        "AF_AFR": ("Allele Count African/African American",
                   "Allele Number African/African American"),
        "AF_AMR": ("Allele Count Admixed American",
                   "Allele Number Admixed American"),
        "AF_EAS": ("Allele Count East Asian",
                   "Allele Number East Asian"),
        "AF_EUR": ("Allele Count European (non-Finnish)",
                   "Allele Number European (non-Finnish)"),
        "AF_SAS": ("Allele Count South Asian",
                   "Allele Number South Asian"),
    }
    for af_name, (count_col, num_col) in pops.items():
        df[count_col] = pd.to_numeric(df[count_col], errors="coerce")
        df[num_col]   = pd.to_numeric(df[num_col],   errors="coerce")
        df[af_name]   = df[count_col] / df[num_col]

    df["AF_popmax"] = df[list(pops.keys())].max(axis=1)

    #1. filter only rsID
    df = df[df["rsIDs"].str.startswith("rs", na=False)]
    #(f"After rsID filter: {len(df)}")

    #2. filter no pLoF (no overlap with Group A)
    pLoF_types = ["stop_gained", "splice_donor", "splice_acceptor",
           "frameshift_variant", "stop_lost", "start_lost"]
    df = df[~df["VEP Annotation"].isin(pLoF_types)]
    #print(f"After pLoF: {len(df)}")

    #3. filter AF > 1% AND < 85% (otherwise no variants was found)
    frequency_filter = (
        ((df["Allele Frequency"] > 0.01) | (df["AF_popmax"] > 0.01)) &
        (df["Allele Frequency"] < 0.85) &
        (df["AF_popmax"] < 0.85)
    )
    df = df[frequency_filter]
    #print(f"After AF filter: {len(df)}")

    #3. filter no pathogenic
    exclude = ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]
    df = df[~df["ClinVar Germline Classification"].isin(exclude)]
    #print(f"After ClinVar filter: {len(df)}")

    #4. exclude V370A (Group C)
    df = df[df["rsIDs"] != "rs3827760"]

    #5. filter coding variants only
    coding = ["missense_variant", "synonymous_variant", "splice_region_variant"]
    df = df[df["VEP Annotation"].isin(coding)]
    #print(f"After coding: {len(df)}")

    # nach allen Filtern, vor der manuellen final_rsids-Auswahl:
    print(df.sort_values("AF_popmax", ascending=False)[["rsIDs", "Protein Consequence", "AF_popmax"]].to_string(
        index=False))

    #then manually select final 6 variants
    final_rsids = [
        "rs61761321",   # M107V  8.4% EAS (missense)
        "rs114808659",  # G239R  6.8% AFR (missense)
        "rs79648056",   # S274S  6.8% AFR (synonymous)
        "rs79798733",   # V246M  3.3% AFR (missense)
        "rs3749099",    # P290P  3.0% EAS (synonymous)
        "rs146567337",  # S380R  2.4% EAS (missense)
    ]

    df_b = df[df["rsIDs"].isin(final_rsids)].copy()

    #manual addition: rs146567337 not in gnomAD_B_all.csv
    if "rs146567337" not in df_b["rsIDs"].values:
        manual = pd.DataFrame([{
            "rsIDs":                          "rs146567337",
            "Protein Consequence":            "p.Ser380Arg",
            "VEP Annotation":                 "missense_variant",
            "Allele Frequency":               0.0008042,
            "AF_EAS":                         0.02403,
            "AF_AFR":                         0.00008,
            "AF_EUR":                         0.000006,
            "AF_SAS":                         0.00047,
            "AF_AMR":                         0.00007,
            "AF_popmax":                      0.02403,
        }])
        df_b = pd.concat([df_b, manual], ignore_index=True)

    df_b["Group"] = "B"
    df_b = df_b.sort_values("AF_popmax", ascending=False).reset_index(drop=True)

    #check for missing variants
    found = set(df_b["rsIDs"].values)
    missing = set(final_rsids) - found
    if missing:
        print(f"Check if missing from file: {missing}")

    df_b.to_csv(OUT_B, index=False)
    print(f"\nGroup B: {len(df_b)} variants saved as {OUT_B}")
    print(df_b[["rsIDs", "Protein Consequence", "AF_popmax"]].to_string(index=False))
    return df_b


#GROUP C: loaded from LD analysis output (with bash)
def load_group_c(path):
    df = pd.read_csv(path)

    #keep only 3 variants: V370A + top 2 LD partners
    keep = ["rs3827760", "rs4676213", "rs72627476"]
    df = df[df["rsID"].isin(keep)].copy()
    df["Group"] = "C"

    print(f"Group C: {len(df)} variants")
    print(df[["rsID", "Protein", "Consequence", "r2_V370A"]].to_string(index=False))
    return df


#VARIANT SET TABLE: Merge A + B + C
def build_variantset(a, b, c):
    a_clean = pd.DataFrame({
        "rsID":        a["rsID"],
        "Protein":     a["Protein_Change"],
        "Consequence": a["Consequence"],
        "ClinVar":     a["ClinVar_Class"],
        "CHR":         a["CHR"],
        "POS_GRCh38":  a["POS_GRCh38"],
        "POS_GRCh37":  pd.NA,
        "AF_EAS":      pd.NA, "AF_AFR": pd.NA, "AF_EUR": pd.NA,
        "AF_SAS":      pd.NA, "AF_AMR": pd.NA, "AF_popmax": pd.NA,
        "r2_V370A":    pd.NA, "Note":   pd.NA, "Group":     "A",
    })

    b_clean = pd.DataFrame({
        "rsID":        b["rsIDs"],
        "Protein":     b["Protein Consequence"],
        "Consequence": b["VEP Annotation"],
        "ClinVar":     b["ClinVar Germline Classification"],
        "CHR":         b.get("Chromosome", 2),
        "POS_GRCh38":  b.get("Position", pd.NA),
        "POS_GRCh37":  pd.NA,
        "AF_EAS":      b["AF_EAS"], "AF_AFR": b["AF_AFR"], "AF_EUR": b["AF_EUR"],
        "AF_SAS":      b["AF_SAS"], "AF_AMR": b["AF_AMR"],
        "AF_popmax":   b["AF_popmax"],
        "r2_V370A":    pd.NA, "Note": pd.NA, "Group": "B",
    })

    c_clean = pd.DataFrame({
        "rsID":        c["rsID"],
        "Protein":     c["Protein"],
        "Consequence": c["Consequence"],
        "ClinVar":     pd.NA,
        "CHR":         2,
        "POS_GRCh38":  pd.NA,
        "POS_GRCh37":  c["POS_GRCh37"],
        "AF_EAS":      c.get("AF_EAS", pd.NA),
        "AF_AFR":      c.get("AF_AFR", pd.NA),
        "AF_EUR":      c.get("AF_EUR", pd.NA),
        "AF_SAS":      pd.NA, "AF_AMR": pd.NA, "AF_popmax": pd.NA,
        "r2_V370A":    c["r2_V370A"], "Note": c["Note"], "Group": "C",
    })

    variantset = pd.concat([a_clean, b_clean, c_clean], ignore_index=True)

    cons_map = {
        "nonsense":                 "stop_gained",
        "missense variant":         "missense_variant",
        "splice donor variant":     "splice_donor_variant",
        "splice acceptor variant":  "splice_acceptor_variant",
    }
    variantset["Consequence"] = variantset["Consequence"].replace(cons_map)

    print(f"After merge: {len(variantset)} variants")
    print(variantset["Consequence"].value_counts().to_string())
    return variantset

if __name__ == "__main__":

    a = build_group_a(CLINVAR_FILE)
    b = build_group_b(GNOMAD_FILE)
    c = load_group_c(GROUPC_FILE)

    variantset = build_variantset(a, b, c)

    variantset.to_csv(OUT_VARIANTSET, index=False)
    print(f"\nSaved as {OUT_VARIANTSET}")