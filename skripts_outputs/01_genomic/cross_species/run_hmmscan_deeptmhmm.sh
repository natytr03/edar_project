#!/bin/bash

set -euo pipefail

DATASET_DIR=~/edar/dataset
WORKDIR=~/edar/domain
INPUT_FASTA="${DATASET_DIR}/species_list_edar_5.fasta"
HUMAN_REF="homo_sapiens"

PFAM_URL="http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

#check if file is correct and create a directionary for output if needed
mkdir -p "$WORKDIR"
cd "$WORKDIR"

log() { echo -e "\n[pipeline] $*"; }

if [ ! -f "$INPUT_FASTA" ]; then
    echo "ERROR: FASTA not found at $INPUT_FASTA" >&2
    exit 1
fi


#1. Pfam-A setup + hmmscan (for DEath Domain and TRNF)
log "STAGE 1: Pfam-A + hmmscan"

if [ ! -f "Pfam-A.hmm" ]; then
    log "CHECK: Downloading Pfam-A.hmm.gz"
    wget -c "$PFAM_URL"
    gunzip -k Pfam-A.hmm.gz
else
    log "Pfam-A.hmm already present. SKIP download."
fi

#if Pfam-A.hmm not there then do hmmpress
if [ ! -f "Pfam-A.hmm.h3m" ]; then
    log "Pfam-A.hmm not yet hmmpress prepared. RUN hmmpress."
    hmmpress Pfam-A.hmm
fi


#check if results already there, if not it searches input sequence against Pfam-A database and save domain matches
DOMTBLOUT="hmmscan_deathdomain.domtblout"
if [ ! -s "$DOMTBLOUT" ]; then
    log "RUN hmmscan against full Pfam-A."
    hmmscan --domtblout "$DOMTBLOUT" \
            Pfam-A.hmm \
            "$INPUT_FASTA" \
            > hmmscan_deathdomain.stdout.log
else
    log "$DOMTBLOUT already exists. SKIP hmmscan."
fi


#2. DeepTMHMM via BioLib  (for signal peptide, ECD, and TM helix)
log "STAGE 2: DeepTMHMM (BioLib)"

TOPOLOGIES_FILE="predicted_topologies.3line"

#run DeepTMHMM if results not there
if [ ! -s "$TOPOLOGIES_FILE" ]; then
    log "ENSURE pybiolib is up to date."
    #make sure to update BioLib
    pip install --upgrade pybiolib --quiet
    log "RUN DeepTMHMM via BioLib cloud (sequences are uploaded to BioLib)."
    (
        cd "$DATASET_DIR"
        biolib run DTU/DeepTMHMM --fasta "$(basename "$INPUT_FASTA")"
    )

    #look for output, only the newest ones
    FOUND_TOPO=$(find "$DATASET_DIR" -name "predicted_topologies.3line" -newer "$INPUT_FASTA" | head -1 || true)
    #if not check again
    if [ -z "$FOUND_TOPO" ]; then
        echo "ERROR: DeepTMHMM run did not produce predicted_topologies.3line. CHECK BioLib output above." >&2
        exit 1
    fi
    #if there is output file move ut to working direction
    mv "$FOUND_TOPO" "$WORKDIR/$TOPOLOGIES_FILE"
    log "Moved DeepTMHMM output to $WORKDIR/$TOPOLOGIES_FILE"
else
    log "$TOPOLOGIES_FILE already exists. SKIP DeepTMHMM run."
fi


#parse output for plot
BOUNDARIES_CSV="domain_boundaries_deeptmhmm.csv"

log "Cross-species domain architecture parsing and plotting"
    python3 parse_and_plot_full_domain_architecture_v2.py \
        --domtblout hmmscan_deathdomain.domtblout \
        --tm-csv tm_domains.csv \
        --out edar_full_domain_architecture_final