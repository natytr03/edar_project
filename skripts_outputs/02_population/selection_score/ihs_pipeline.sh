#!/bin/bash

set -e

#define paths:
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

INPUTDIR="$PROJECT_ROOT/data"
WORKDIR="$PROJECT_ROOT/skripts_outputs/02_population/selection_score"

# iput files:
VCF_ALL="$INPUTDIR/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
EAS_SAMPLES="$INPUTDIR/EAS_samples.txt"
GMAP="$INPUTDIR/plink.chr2.GRCh37.map"

# define region for scan
REGION="2:107000000-112000000"
V370A_POS="109513601"

cd "$WORKDIR" || exit 1

#create eas vcf, as iHS calculation in eas
EAS_VCF="$INPUTDIR/edar_eas_5mb.vcf.gz"

if [ ! -f "$EAS_VCF" ]; then
    tabix -h "$VCF_ALL" "$REGION" \
    | bcftools view \
        -S "$EAS_SAMPLES" \
        --min-ac 1 \
        --min-alleles 2 \
        --max-alleles 2 \
        --type snps \
        -o "$EAS_VCF" \
        -O z

    tabix -p vcf "$EAS_VCF"
    echo "vcf eas is created"
else
    echo "skip, vsf already there"
fi

#short overview check vcf file
N_SAMPLES=$(bcftools query -l "$EAS_VCF" | wc -l)
N_SNPS=$(bcftools stats "$EAS_VCF" | grep "^SN.*SNPs" | cut -f4)

echo "  Samples: $N_SAMPLES | SNPs: $N_SNPS"

#for each point SNP positions from EAS-VCF and add an interpolated genetic position in cM
#creating a .gmap file needed for selscan
if [ ! -f edar_5mb_full.gmap ]; then

    bcftools query -f '%POS\n' "$INPUTDIR/edar_eas_5mb.vcf.gz" > vcf_5mb_positions.txt

    python3 << EOF
import pandas as pd
import numpy as np

map_df = pd.read_csv("$GMAP",sep="\t", header=None, names=["chr", "id", "cM", "bp"])

vcf_pos = pd.read_csv("$WORKDIR/vcf_5mb_positions.txt", header=None, names=["bp"])

vcf_pos["cM"] = np.interp(vcf_pos["bp"],
                          map_df["bp"],
                          map_df["cM"])

# selscan Format: CHR  ID  cM  bp
vcf_pos["chr"] = 2
vcf_pos["id"] = "."

vcf_pos[["chr", "id", "cM", "bp"]].to_csv(
    "$WORKDIR/edar_5mb_full.gmap",
    sep="\t", index=False, header=False
)

print(f"Map created")
EOF
    echo ".gmap file is created"
else
    echo "skip, .gmap file already there"
fi

#iHS calcuation
GMAP_OUT="$WORKDIR/edar_5mb_full.gmap"

selscan \
    --ihs \
    --vcf "$EAS_VCF" \
    --map "$GMAP_OUT" \
    --out "$IHS_PREFIX" \
    --threads 4

N_LOCI=$(wc -l < $WORKDIR/edar_eas_5mb_ihs.ihs.out)

#for comparison iHS results need to be normalized
norm --ihs \
    --files $WORKDIR/edar_eas_5mb_ihs.ihs.out

#create txt file for output results
RESULTS_TXT="$WORKDIR/edar_eas_5mb_ihs_results.txt"

{
    #check if V370A is in raw iHS output
    echo "V370A check in raw iHS:"
    grep -w "$V370A_POS" "$WORKDIR/edar_eas_5mb_ihs.ihs.out" \
        || echo "V370A is not in output"

    #normalized iHS file
    NORM_FILE="$WORKDIR/edar_eas_5mb_ihs.ihs.out.100bins.norm"

    echo ""
    echo "V370A (rs3827760, pos $V370A_POS):"

    grep -w "$V370A_POS" "$NORM_FILE" \
        | awk '{printf "  pos=%s  freq=%.3f  iHS_raw=%.3f  iHS_norm=%.3f  flag=%s\n", $2, $3, $6, $7, $8}'

    echo ""
    echo "Top 10 signals in this region"
    echo "  pos          freq    iHS_raw  iHS_norm  flag"

    # Sort by absolute iHS_norm
    awk '{abs = ($7 < 0 ? -$7 : $7)
        print abs, $0}' "$NORM_FILE" \
        | sort -k1,1nr \
        | head -10 \
        | cut -d' ' -f2- \
        | awk '{printf "  %-12s %.3f   %.3f    %.3f     %s\n",
                $2, $3, $6, $7, $8}'

    echo ""
    echo "Significant loci (|iHS_norm| > 2)"

    #count loci with absolute normalized iHS > 2
    N_SIG=$(awk '$7 != "NA" && ($7 > 2 || $7 < -2)' "$NORM_FILE" | wc -l)

    N_TOTAL=$(wc -l < "$NORM_FILE")

    PERCENT=$(echo "scale=1; $N_SIG*100/$N_TOTAL" | bc)

    echo "  $N_SIG von $N_TOTAL Loci (${PERCENT}%)"

    echo ""
    echo "Saved Outputpaths:"
    echo "  Raw iHS:  $WORKDIR/edar_eas_5mb_ihs.ihs.out"
    echo "  Norm iHS: $NORM_FILE"
    echo "  Results:  $RESULTS_TXT"

    echo ""
    echo "Saved."

} | tee "$RESULTS_TXT"