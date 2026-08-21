from Bio import Entrez, SeqIO
import time
import csv
import re
from datetime import datetime

#konfiguration
Entrez.email   = "nltran468@gmail.com"
OUTPUT_FASTA   = "species_list_edar_cds.fasta"
SLEEP_SEC      = 0.4

SEQUENCES = [
    # (FASTA-Label,                                          Accession,           Gruppe,           Sonderfall)
    ("homo_sapiens_edar",                                    "NM_022336.4",       "Primates",        None),
    ("pan_troglodytes_edar",                                 "XM_001139583.7",    "Primates",        None),
    ("gorilla_gorilla_edar",                                 "XM_055378332.2",    "Primates",        None),
    ("ursus_arctos_edar",                                    "XM_026516135.4",    "Hairy_mammals",   None),
    ("bos_mutus_edar",                                       "XM_070380157.1",    "Hairy_mammals",   None),
    ("manis_javanica_edar",                                  "XM_037015587.2",    "Hairless_mammals",None),
    ("neophocaena_asiaeorientalis_edar",                     "XM_024757981.1",    "Hairless_mammals",None),
    ("loxodonta_africana_edar",                              "XM_023552651.2",    "Hairless_mammals",None),
    ("trichechus_manatus_latirostris_edar",                  "XM_023732242.1",    "Hairless_mammals",None),
    ("sus_scrofa_edar",                                      "XM_013995872.2",    "Hairy_mammals",   None),
    ("heterocephalus_glaber_edar",                           "XM_004844733.2",    "Hairless_mammals",None),
    ("struthio_camelus_edar",                                "XM_009665487.2",    "Birds",           None),
    ("dromaius_novaehollandiae_edar",                        "XM_026097008.2",    "Birds",           None),
    ("anser_cygnoides_edar",                                 "XM_066990475.1",    "Birds",           None),
    ("columba_livia_edar",                                   "XM_065059922.1",    "Birds",           None),
    ("alligator_mississippiensis_edar",                      "XM_019485001.2",    "Reptiles",        None),
    ("chrysemys_picta_bellii_edar",                          "XM_024110836.3",    "Reptiles",        None),
    ("python_bivittatus_edar",                               "XM_025166306.1",    "Reptiles",        None),
    ("pogona_vitticeps_edar",                                "XM_078391467.1",    "Reptiles",        None),
    ("danio_rerio_edar",                                     "NM_001115064.3",    "Fish_outgroup",   None),
    ("cyprinus_carpio_edar",                                 "XM_042755187.1",    "Fish_outgroup",   None),
    ("oncorhynchus_clarkii_lewisi_edar",                     "XM_071122972.1",    "Fish_outgroup",   None),
    ("scyliorhinus_canicula_edar",                           "XM_038819029.1",    "Fish_outgroup",   None),
    ("archocentrus_centrarchus_edar",                        "XM_030757580.1",    "Fish_outgroup",   None),
]



def get_db_for_accession(accession):

    if  accession.startswith(("NM_", "XM_", "NR_")):
        return "nucleotide"
    elif re.match(r"^[A-Z][0-9]{5}", accession):  #UniProt format
        return "protein"
    else:
        return "nucleotide"


def process_entry(label, accession, group):

    print(f"\n  [{group}] {label}")

    #Accession fetchen
    db = get_db_for_accession(accession)
    print(f"    Accession: {accession} → DB: {db}")

    try:
            #if nucleotide accession
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            time.sleep(SLEEP_SEC)

            cds_seq = None
            for feature in record.features:
                #find cds feature and save for fasta
                if feature.type == "CDS":
                    cds_seq = str(feature.extract(record.seq))
                    break #only take first one

            if not cds_seq:
                print(f"No CDS found. {record.id}")
                return {"label": label, "group": group, "accession": accession,
                        "status": "NO_CDS", "cds_length": 0, "cds_seq": ""}

    except Exception as e:
        print(f"ERROR: {e}")
        return {"label": label, "group": group, "accession": accession,
                "status": f"EXCEPTION: {e}", "cds_length": 0, "cds_seq": ""}

    #check length
    cds_length = len(cds_seq)

    if cds_length % 3 != 0:
        length_flag= "need to check: CDS with unusual length"
    else:
        length_flag = "CDS with expected length"

    print(f" {length_flag} {cds_length} bp")
    return {"label": label, "group": group, "accession": accession,
            "status": "OK", "cds_length": cds_length, "cds_seq": cds_seq}


def main():

    all_results = []

    for label, accession, group, special in SEQUENCES:
        result = process_entry(label, accession, group)
        all_results.append(result)

    ok     = [r for r in all_results if r["status"] == "OK"]
    #save as fasta format
    with open(OUTPUT_FASTA, "w") as f:
        for r in ok:
            f.write(f">{r['label']}\n{r['cds_seq']}\n\n")
    print(f"\n  OUTPUT_FASTA: {OUTPUT_FASTA}  ({len(ok)} Sequences)")

if __name__ == "__main__":
    main()