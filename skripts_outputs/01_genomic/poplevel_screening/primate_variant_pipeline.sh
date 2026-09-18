#!/usr/bin/env bash

set -euo pipefail

WORKDIR="${1:-edar_pipeline_run}"
mkdir -p "$WORKDIR"; cd "$WORKDIR"

HUMAN_PROTEIN_ID="ENSP00000258443"
GAPDH_PROTEIN_ID="ENSP00000229239"
MACAQUE_PROTEIN_ID="ENSMMUP00000036157"
TARGET_AA_POS=370
RSID="rs3827760"
GREAT_APE_CHROM="chr2"
FLANK_WINDOW=2000
LIFTOVER_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

declare -A VCF_URLS=(
 [Gorilla]="https://app05b.phaidra.org/pfsa/o_2066302/merged_segregating/Gorilla/Gorilla_wild_filtered/chr2.filteranno.vcf.gz"
 [Pan]="https://app05b.phaidra.org/pfsa/o_2066302/merged_segregating/Pan/Pan_wild_filtered/chr2.filteranno.vcf.gz"
 [Pongo]="https://app05b.phaidra.org/pfsa/o_2066302/merged_segregating/Pongo/Pongo_wild_filtered/chr2.filteranno.vcf.gz"
)
declare -A GIBBON_URLS=(
 [Nomascus_siki]="https://download.cncb.ac.cn/gwh/Animals/Nomascus_siki_NSI_Genome_GWHFIUK00000000.1/GWHFIUK00000000.1.genome.fasta.gz"
 [Hoolock_leuconedys]="https://download.cncb.ac.cn/gwh/Animals/Hoolock_leuconedys_HLE_Genome_GWHFIUM00000000.1/GWHFIUM00000000.1.genome.fasta.gz"
 [Symphalangus_syndactylus]="https://download.cncb.ac.cn/gwh/Animals/Symphalangus_syndactylus_SSY_Genome_GWHFIUJ00000000.1/GWHFIUJ00000000.1.genome.fasta.gz"
 [Hylobates_pileatus]="https://download.cncb.ac.cn/gwh/Animals/Hylobates_pileatus_HPI_Genome_GWHFIUL00000000.1/GWHFIUL00000000.1.genome.fasta.gz"
)

die(){ echo "ERROR: $*" >&2; exit 1; }
warn(){ echo "WARNING: $*" >&2; }

download(){
  local url=$1 dest=$2 min=${3:-1}
  [[ -s "$dest" ]] && [[ $(stat -c%s "$dest" 2>/dev/null || stat -f%z "$dest") -ge $min ]] && {
    [[ "$dest" != *.gz ]] || gzip -t "$dest"
    return
  }
  rm -f "$dest"
  curl -fL --retry 3 --retry-delay 5 --connect-timeout 30 -o "$dest" "$url" || {
    rm -f "$dest"; return 1;
  }
  [[ "$dest" != *.gz ]] || gzip -t "$dest" || { rm -f "$dest"; return 1; }
}

get_protein(){
  local id=$1 out=$2
  [[ -s "$out" ]] && return
  curl -fLsS --retry 3 --retry-delay 5 --connect-timeout 30 --max-time 120 \
    "https://rest.ensembl.org/sequence/id/$id" \
    -H "Content-Type: text/x-fasta" -o "$out" ||
    die "Could not download $id"
  grep -q '^>' "$out" || die "Invalid FASTA: $out"
}

human_proteins(){
  get_protein "$HUMAN_PROTEIN_ID" human_edar.fasta
  get_protein "$GAPDH_PROTEIN_ID" human_gapdh.fasta
}


human_coordinates(){
  download "$LIFTOVER_URL" hg19ToHg38.over.chain.gz 1000 ||
    die "Could not obtain liftover chain"

  python3 - "$RSID" > human_coords.json <<'PY'
import sys
import json
import requests

rsid = sys.argv[1]
def get_mapping(base_url):
    url = f"{base_url}/variation/human/{rsid}"
    r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=60)
    r.raise_for_status()
    data = r.json()
    standard = [str(i) for i in range(1, 23)] + ["X", "Y"]
    return [m for m in data.get("mappings", [])
    if str(m.get("seq_region_name")) in standard]

grch37 = get_mapping("https://grch37.rest.ensembl.org")
grch38 = get_mapping("https://rest.ensembl.org")

if not grch37:
    raise RuntimeError("No GRCh37 mapping found")

if not grch38:
    raise RuntimeError("No GRCh38 mapping found")

result = {
    "rsid": rsid,
    "chrom": str(grch37[0]["seq_region_name"]),
    "grch37_pos": int(grch37[0]["start"]),
    "grch38_pos": int(grch38[0]["start"]),
    "allele_string": grch38[0].get("allele_string")
}

print(json.dumps(result, indent=2))
PY

  python3 <<'PY'
import json
from pyliftover import LiftOver

with open("human_coords.json") as f:
    coords = json.load(f)

lo = LiftOver("hg19ToHg38.over.chain.gz")
chrom = str(coords["chrom"])
pos0 = int(coords["grch37_pos"]) - 1

if not chrom.startswith("chr"):
    chrom = "chr" + chrom

result = lo.convert_coordinate(chrom, pos0)

if not result:
    raise SystemExit("Liftover returned no result")

lifted_chr = result[0][0]
lifted_pos = result[0][1] + 1

expected_chr = "chr" + str(coords["chrom"])
expected_pos = int(coords["grch38_pos"])

match = (
    lifted_chr == expected_chr
    and lifted_pos == expected_pos
)

print(
    f"Liftover cross-check: "
    f"Ensembl GRCh38={expected_chr}:{expected_pos}, "
    f"liftover={lifted_chr}:{lifted_pos} "
    f"{'OK' if match else 'MISMATCH'}"
)

if not match:
    raise SystemExit("Liftover mismatch")
PY

  GRCH38_POS=$(python3 -c "import json;print(json.load(open('human_coords.json'))['grch38_pos'])")
  HUMAN_CHROM="chr$(python3 -c "import json;print(json.load(open('human_coords.json'))['chrom'])")"
  export GRCH38_POS HUMAN_CHROM

  python3 <<'PY'
import json
d=json.load(open("human_coords.json"))
with open("human_target_coordinate.txt","w") as f:
    f.write(f"RSID={d['rsid']}\n")
    f.write(f"GRCh37=chr{d['chrom']}:{d['grch37_pos']}\n")
    f.write(f"GRCh38=chr{d['chrom']}:{d['grch38_pos']}\n")
    f.write(f"Alleles={d.get('allele_string')}\n")
PY
}

great_apes(){
  local region="$GREAT_APE_CHROM:$GRCH38_POS-$GRCH38_POS"
  local ws=$((GRCH38_POS-FLANK_WINDOW))
  local we=$((GRCH38_POS+FLANK_WINDOW))
  ((ws<1)) && ws=1
  local window="$GREAT_APE_CHROM:$ws-$we"
  echo "clade,target_records,window_records,status" > great_apes_summary.csv

  for clade in Gorilla Pan Pongo; do
    local vcf="${clade}_chr2.vcf.gz"
    download "${VCF_URLS[$clade]}" "$vcf" 10000000 || {
      echo "$clade,NA,NA,DOWNLOAD_FAILED" >> great_apes_summary.csv; continue;
    }
    bcftools index -f -t "$vcf" || {
      echo "$clade,NA,NA,INDEX_FAILED" >> great_apes_summary.csv; continue;
    }

    local target="${clade}_target.txt" flank="${clade}_window.txt"
    bcftools view -r "$region" "$vcf" >"$target" 2>/dev/null || {
      echo "$clade,NA,NA,QUERY_FAILED" >> great_apes_summary.csv; continue;
    }
    bcftools view -r "$window" "$vcf" >"$flank" 2>/dev/null || {
      echo "$clade,NA,NA,WINDOW_QUERY_FAILED" >> great_apes_summary.csv; continue;
    }

    local nt nw
    nt=$(grep -vc '^#' "$target" || true)
    nw=$(grep -vc '^#' "$flank" || true)
    if ((nt>0)); then
      echo "$clade,$nt,$nw,VARIANT_FOUND" >> great_apes_summary.csv
    else
      echo "$clade,0,$nw,NO_VARIANT" >> great_apes_summary.csv
    fi
  done
}

macaque(){
  python3 - "$MACAQUE_PROTEIN_ID" "$TARGET_AA_POS" > macaque_coords.txt <<'PY'
import sys
import requests
from Bio import Align
from Bio.Align import substitution_matrices

HUMAN_PROTEIN_ID = "ENSP00000258443"
macaque_id = sys.argv[1]
target_pos = int(sys.argv[2])

def get_seq(protein_id):
    url = f"https://rest.ensembl.org/sequence/id/{protein_id}"
    r = requests.get(url, params={"type": "protein", "content-type": "text/x-fasta"}, timeout=60)
    r.raise_for_status()
    lines = r.text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError(f"Unexpected FASTA response for {protein_id}")
    return "".join(lines[1:])

human_seq = get_seq(HUMAN_PROTEIN_ID)
macaque_seq = get_seq(macaque_id)

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.mode = "global"

alignment = aligner.align(human_seq, macaque_seq)[0]
q = str(alignment[0])
s = str(alignment[1])

q_idx = 0
s_idx = 0
macaque_pos = None

for qc, sc in zip(q, s):
    if qc != "-":
        q_idx += 1
    if sc != "-":
        s_idx += 1
    if q_idx == target_pos:
        if sc != "-":
            macaque_pos = s_idx
        break

if macaque_pos is None:
    raise RuntimeError(f"Human EDAR AA {target_pos} does not map to a macaque amino-acid residue")

url = f"https://rest.ensembl.org/map/translation/{macaque_id}/{macaque_pos}..{macaque_pos}"
r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=60)
r.raise_for_status()
data = r.json()

if not data.get("mappings"):
    raise RuntimeError("Ensembl returned no macaque genomic mapping")

m = data["mappings"][0]

print(f"human_protein={HUMAN_PROTEIN_ID}")
print(f"macaque_protein={macaque_id}")
print(f"human_target_aa={target_pos}")
print(f"macaque_protein_pos={macaque_pos}")
print(f"assembly={m['assembly_name']}")
print(f"chrom={m['seq_region_name']}")
print(f"start={m['start']}")
print(f"end={m['end']}")
PY
}

gibbon_result(){
  python3 - "$1" "$TARGET_AA_POS" <<'PY'
import sys
path,pos=sys.argv[1],int(sys.argv[2]); calls=set()
for line in open(path):
    f=line.rstrip().split("\t")
    if len(f)<10: continue
    try: qs,qe=map(int,f[4:6])
    except ValueError: continue
    if not qs<=pos<=qe: continue
    qi=qs-1
    for qc,sc in zip(f[8],f[9]):
        if qc!="-": qi+=1
        if qi==pos:
            if sc!="-": calls.add(sc)
            break
if not calls: print("NO_COVERAGE")
elif len(calls)>1: print("AMBIGUOUS:"+",".join(sorted(calls)))
else: print(next(iter(calls)))
PY
}

gibbons(){
  echo "species,result,gapdh_control" > gibbons_summary.csv
  human_proteins

  for species in Nomascus_siki Hoolock_leuconedys Symphalangus_syndactylus Hylobates_pileatus; do
    local gz="${species}_genome.fasta.gz" fasta="${species}_genome.fasta" db="${species}_db"
    download "${GIBBON_URLS[$species]}" "$gz" 50000000 || {
      echo "$species,DOWNLOAD_FAILED,NOT_TESTED" >> gibbons_summary.csv; continue;
    }
    [[ -s "$fasta" ]] || gunzip -f "$gz" || {
      echo "$species,GZIP_FAILED,NOT_TESTED" >> gibbons_summary.csv; continue;
    }
    [[ -s "$fasta" ]] || { echo "$species,EMPTY_FASTA,NOT_TESTED" >> gibbons_summary.csv; continue; }

    [[ -f "${db}.nin" || -f "${db}.00.nin" ]] ||
      makeblastdb -in "$fasta" -dbtype nucl -parse_seqids -out "$db" || {
        echo "$species,BLASTDB_FAILED,NOT_TESTED" >> gibbons_summary.csv; continue;
      }

    local hits="${species}_edar_hits.tsv"
    tblastn -query human_edar.fasta -db "$db" \
      -outfmt "6 qseqid sseqid pident length qstart qend sstart send qseq sseq" \
      -max_target_seqs 3 -max_hsps 3 -evalue 1e-5 >"$hits" || {
        echo "$species,TBLASTN_FAILED,NOT_TESTED" >> gibbons_summary.csv; continue;
      }

    local result gapdh=NOT_NEEDED
    result=$(gibbon_result "$hits")

    if [[ "$result" == NO_COVERAGE ]]; then
      local gf="${species}_gapdh_control.txt"
      if tblastn -query human_gapdh.fasta -db "$db" -outfmt 0 \
          -max_target_seqs 1 -evalue 1e-5 >"$gf" &&
         grep -q "Sequences producing significant alignments" "$gf"; then
        gapdh=SUCCESS
      else
        gapdh=FAILED
      fi
    fi
    echo "$species,$result,$gapdh" >> gibbons_summary.csv
  done
}

summary(){
  {
    echo "EDAR primate variant pipeline"
    echo "RSID: $RSID"
    echo "Human protein: $HUMAN_PROTEIN_ID"
    echo "Target amino acid: $TARGET_AA_POS"
    echo "Human target: $HUMAN_CHROM:$GRCH38_POS"
    echo
    echo "Great Apes"; cat great_apes_summary.csv
    echo; echo "Macaque"; cat macaque_coords.txt
    echo; echo "Gibbons"; cat gibbons_summary.csv
    echo; echo "Files generated:"
    find . -maxdepth 1 -type f -printf '%f\n' | sort
  } > pipeline_summary.txt
}

main(){
  human_proteins
  human_coordinates
  great_apes
  macaque
  gibbons
  summary
  echo "Pipeline completed"
}
main "$@"

