import pandas as pd
import requests
import time


#candidates show clear gap, where two have high r² (0.886–0.897), while the others are lower (≤0.79)
#full candidate list saved in EDAR_Gruppe_C_all_candidates.csv
R2_THRESHOLD = 0.88

#get the rsIDs from position
def get_rsid(pos, chrom=2, server="https://grch37.rest.ensembl.org"):
    url = f"/overlap/region/human/{chrom}:{pos}-{pos}?feature=variation;content-type=application/json"
    req = requests.get(server + url)
    if req.ok:
        data = req.json()
        rsids = [v["id"] for v in data if v["id"].startswith("rs")]
        return rsids[0] if rsids else "."
    return "."

#get the corresponding consequence of the variant
def get_consequence(rsid, server="https://grch37.rest.ensembl.org"):
    if rsid == ".":
        return "unknown"
    url = f"/vep/human/id/{rsid}?content-type=application/json"
    req = requests.get(server + url)
    if req.ok:
        data = req.json()
        if data and "most_severe_consequence" in data[0]:
            return data[0]["most_severe_consequence"]
    return "unknown"


df = pd.read_csv(r"EDAR_Gruppe_C.csv")

#run through every variant to get the information needed
is_ld_partner = df["rsID"] == "."
for idx in df[is_ld_partner].index:
    pos = df.loc[idx, "POS_GRCh37"]
    rsid = get_rsid(pos)
    df.loc[idx, "rsID"] = rsid
    print(f"chr2:{pos} -> {rsid}")
    time.sleep(0.5)

    consequence = get_consequence(rsid)
    df.loc[idx, "Consequence"] = consequence
    print(f"  {rsid} -> {consequence}")
    time.sleep(0.5)


df.to_csv("EDAR_Gruppe_C_all_candidates.csv", index=False)
print(f"\nGespeichert: EDAR_Gruppe_C_all_candidates.csv")

#save the final list of group c
gruppe_c_final = df[(df["rsID"] == "rs3827760") | (df["r2_V370A"] > R2_THRESHOLD)].copy()
gruppe_c_final = gruppe_c_final.sort_values("r2_V370A", ascending=False).reset_index(drop=True)

gruppe_c_final.to_csv("EDAR_Gruppe_C_new.csv", index=False)