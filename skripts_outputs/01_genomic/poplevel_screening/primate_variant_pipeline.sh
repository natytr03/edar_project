#!/usr/bin/env bash

set -uo pipefail

WORKDIR="${1:-edar_pipeline_run}"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

HUMAN_PROTEIN_ID="ENSP00000258443"
GAPDH_PROTEIN_ID="ENSP00000229239"
MACAQUE_PROTEIN_ID="ENSMMUP00000036157"
TARGET_AA_POS=370
RSID="rs3827760"
GREAT_APE_CHROM="chr2"
FLANK_WINDOW=2000

log()  { echo "[$(date +%H:%M:%S)] $*"; }
warn() { echo "[$(date +%H:%M:%S)] WARNING: $*" >&2; }

#download needed files with clean resume otherwise file became corrupt
download() {
  local url="$1" dest="$2" min_bytes="${3:-1024}"

  if [[ -f "$dest" ]] && [[ $(stat -c%s "$dest" 2>/dev/null || stat -f%z "$dest") -ge "$min_bytes" ]]; then
    log "files already present ($dest)"
    return 0
  fi

  curl -C - -o "$dest" "$url"

  if [[ "$dest" == *.gz ]]; then
    if ! file "$dest" | grep -q "gzip compressed data"; then
      warn "$dest is not valid gzip data. NEED RESTART"
      rm -f "$dest"
      curl -o "$dest" "$url"
    fi
  fi
}


#check first human for correct mapping (in GRCh37/38, with liftover cross-check)
human_coordinates_check() {
    python3 - "$RSID" > stage1_human_coords.json <<'PYEOF'
import sys, json, requests

rsid = sys.argv[1]

def get_mapping(base_url):
    url = f"{base_url}/variation/human/{rsid}"
    r = requests.get(url, headers={"Content-Type": "application/json"})
    r.raise_for_status()
    standard = [str(i) for i in range(1, 23)] + ["X", "Y"]
    return [m for m in r.json()["mappings"] if m["seq_region_name"] in standard]

grch37 = get_mapping("https://grch37.rest.ensembl.org")
grch38 = get_mapping("https://rest.ensembl.org")

result = {
    "chrom": grch37[0]["seq_region_name"],
    "grch37_pos": grch37[0]["start"],
    "grch38_pos": grch38[0]["start"],
    "allele_string": grch38[0]["allele_string"],
}
print(json.dumps(result, indent=2))
PYEOF

  #liftover cross-check
  python3 - <<'PYEOF'
import json
from pyliftover import LiftOver

with open("stage1_human_coords.json") as f:
    coords = json.load(f)

lo = LiftOver("hg19", "hg38")
result = lo.convert_coordinate(f"chr{coords['chrom']}", coords["grch37_pos"] - 1)
lifted_pos38 = result[0][1] + 1 if result else None

match = (lifted_pos38 == coords["grch38_pos"])
print(f"Liftover cross-check: Ensembl GRCh38={coords['grch38_pos']}, "
      f"liftover={lifted_pos38} -> {'OK' if match else 'MISMATCH - INVESTIGATE'}")
PYEOF

  GRCH38_POS=$(python3 -c "import json; print(json.load(open('stage1_human_coords.json'))['grch38_pos'])")
  log "Human target coordinate: ${GREAT_APE_CHROM}:${GRCH38_POS} (GRCh38)"
  export GRCH38_POS

}


#ckeck in greate ape
great_apes_check() {
  local region="${GREAT_APE_CHROM}:${GRCH38_POS}-${GRCH38_POS}"
  local window_region="${GREAT_APE_CHROM}:$((GRCH38_POS - FLANK_WINDOW))-$((GRCH38_POS + FLANK_WINDOW))"

  declare -A VCF_URLS=(
    [Gorilla]="https://app05b.phaidra.org/pfsa/o_2066302/merged_segregating/Gorilla/Gorilla_wild_filtered/chr2.filteranno.vcf.gz"
    [Pan]="https://app05b.phaidra.org/pfsa/o_2066302/merged_segregating/Pan/Pan_wild_filtered/chr2.filteranno.vcf.gz"
    [Pongo]="https://app05b.phaidra.org/pfsa/o_2066302/merged_segregating/Pongo/Pongo_wild_filtered/chr2.filteranno.vcf.gz"
  )

  for clade in "${!VCF_URLS[@]}"; do
    log "$clade"
    local url="${VCF_URLS[$clade]}"
    local target_out window_out source

    #check remote is possible else local download
    target_out=$(bcftools view -r "$region" "$url" 2>/dev/null | grep -v "^#")

    if [[ $? -ne 0 ]]; then
      warn "remote streaming error ($clade). LOCAL DOWNLOAD"
      download "$url" "${clade}_chr2.vcf.gz" 10000000
      download "${url}.csi" "${clade}_chr2.vcf.gz.csi" 1000
      source="${clade}_chr2.vcf.gz"
    else
      source="$url"
    fi

    target_out=$(bcftools view -r "$region" "$source" 2>/dev/null | grep -v "^#")
    window_out=$(bcftools view -r "$window_region" "$source" 2>/dev/null | grep -v "^#")

    local n_target n_window
    n_target=$(echo -n "$target_out" | grep -c . || true)
    n_window=$(echo -n "$window_out" | grep -c . || true)

    echo "${clade},${n_target},${n_window}" >> great_apes_summary.csv

    if [[ "$n_target" -gt 0 ]]; then
      log "  $clade: ALANINE ALLELE FOUND at target position!"
    else
      log "  $clade: no variant at target position (flanking variants in window: $n_window)"
      if [[ "$n_window" -eq 0 ]]; then
        warn "  No flanking variants"
      fi
    fi
  done
}

#check in macaque and mapping on Mul_10 (manual mGAP check in their database as well)
macaque_check() {

  if [[ ! -f human_edar.fasta ]]; then
    curl -s "https://rest.ensembl.org/sequence/id/${HUMAN_PROTEIN_ID}" \
      -H "Content-Type: text/x-fasta" > human_edar.fasta
  fi

  python3 - "$MACAQUE_PROTEIN_ID" "$TARGET_AA_POS" <<'PYEOF'
import sys, requests
from Bio import Align
from Bio.Align import substitution_matrices

macaque_id, target_pos = sys.argv[1], int(sys.argv[2])

def get_seq(protein_id):
    r = requests.get(f"https://rest.ensembl.org/sequence/id/{protein_id}",
                      params={"type": "protein", "content-type": "text/x-fasta"})
    r.raise_for_status()
    return "".join(r.text.strip().splitlines()[1:])

human_seq = get_seq("ENSP00000258443")
macaque_seq = get_seq(macaque_id)

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.mode = "global"
alignment = aligner.align(human_seq, macaque_seq)[0]
q, s = str(alignment[0]), str(alignment[1])

q_idx = s_idx = 0
macaque_pos = None
for qc, sc in zip(q, s):
    if qc != "-": q_idx += 1
    if sc != "-": s_idx += 1
    if q_idx == target_pos:
        macaque_pos = s_idx
        break

r = requests.get(f"https://rest.ensembl.org/map/translation/{macaque_id}/{macaque_pos}..{macaque_pos}",
                  headers={"Content-Type": "application/json"})
r.raise_for_status()
m = r.json()["mappings"][0]
print(f"macaque_protein_pos={macaque_pos}")
print(f"assembly={m['assembly_name']}")
print(f"chrom={m['seq_region_name']}")
print(f"start={m['start']}")
print(f"end={m['end']}")
PYEOF
}

#for gibbons genome tblast was used with human sequence as reference
gibbons_genome_check() {

  if [[ ! -f human_gapdh.fasta ]]; then
    curl -s "https://rest.ensembl.org/sequence/id/${GAPDH_PROTEIN_ID}" \
      -H "Content-Type: text/x-fasta" > human_gapdh.fasta
  fi

  declare -A GIBBON_URLS=(
    [Nomascus_siki]="https://download.cncb.ac.cn/gwh/Animals/Nomascus_siki_NSI_Genome_GWHFIUK00000000.1/GWHFIUK00000000.1.genome.fasta.gz"
    [Hoolock_leuconedys]="https://download.cncb.ac.cn/gwh/Animals/Hoolock_leuconedys_HLE_Genome_GWHFIUM00000000.1/GWHFIUM00000000.1.genome.fasta.gz"
    [Symphalangus_syndactylus]="https://download.cncb.ac.cn/gwh/Animals/Symphalangus_syndactylus_SSY_Genome_GWHFIUJ00000000.1/GWHFIUJ00000000.1.genome.fasta.gz"
    [Hylobates_pileatus]="https://download.cncb.ac.cn/gwh/Animals/Hylobates_pileatus_HPI_Genome_GWHFIUL00000000.1/GWHFIUL00000000.1.genome.fasta.gz"
  )

  for species in "${!GIBBON_URLS[@]}"; do
    log "$species"
    local url="${GIBBON_URLS[$species]}"
    local gz="${species}_genome.fasta.gz"
    local fasta="${species}_genome.fasta"
    local db="${species}_db"

    if [[ ! -f "$fasta" ]]; then
      download "$url" "$gz" 50000000
      gunzip -f "$gz"
    fi

    if [[ ! -f "${db}.nsq" ]] && [[ ! -f "${db}.00.nsq" ]]; then
      makeblastdb -in "$fasta" -dbtype nucl -out "$db"
    fi

    tblastn -query human_edar.fasta -db "$db" \
      -outfmt "6 qseqid sseqid pident length qstart qend sstart send qseq sseq" \
      -max_target_seqs 3 -max_hsps 3 > "${species}_edar_hits.tsv"

    #extract the amino acid at TARGET_AA_POS
    RESULT=$(python3 - "${species}_edar_hits.tsv" "$TARGET_AA_POS" <<'PYEOF'
import sys

path, target_pos = sys.argv[1], int(sys.argv[2])
calls = set()
with open(path) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            continue
        qstart, qend = int(fields[4]), int(fields[5])
        qseq, sseq = fields[8], fields[9]
        if not (qstart <= target_pos <= qend):
            continue
        q_idx = qstart - 1
        for qc, sc in zip(qseq, sseq):
            if qc != "-":
                q_idx += 1
            if q_idx == target_pos:
                calls.add(sc)
                break

if not calls:
    print("NO_COVERAGE")
elif len(calls) > 1:
    print("AMBIGUOUS:" + ",".join(sorted(calls)))
else:
    print(list(calls)[0])
PYEOF
)

    if [[ "$RESULT" == "NO_COVERAGE" ]]; then
      warn "  No BLAST coverage of target position for $species. Cross-check GAPDH positive control."
      tblastn -query human_gapdh.fasta -db "$db" -outfmt 0 -max_target_seqs 1 \
        > "${species}_gapdh_control.txt"
      if grep -q "Sequences producing significant alignments" "${species}_gapdh_control.txt"; then
        log "  Status: unresolved (GAPDH control succeeded)"
      else
        warn "  Status: unresolved (GAPDH control failed)"
      fi
    elif [[ "$RESULT" == AMBIGUOUS:* ]]; then
      warn "  $species: disagreeing amino acids across HSPs ($RESULT)"
    else
      log "  $species: amino acid at target position: $RESULT"
    fi

    echo "${species},${RESULT}" >> gibbons_summary.csv
  done
}


main() {
  : > great_apes_summary.csv
  : > gibbons_summary.csv

  human_coordinates_check
  great_apes_check
  macaque_check | tee macaque_coords.txt
  gibbons_genome_check

  log "Saved."
}

main "$@"



