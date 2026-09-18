#!/usr/bin/env bash

set -euo pipefail

WT_PDB="${1:?Usage: $0 <wt_pdb> <mut_pdb> <position> [chain] [output_name]}"
MUT_PDB="${2:?Missing mutant PDB file}"
POSITION="${3:?Missing target position}"
CHAIN="${4:-A}"
OUT_NAME="${5:-pos${POSITION}}"

#setting default parameters for my analysis
RADII="${RADII:-5 8 10}"
DOMAINS="${DOMAINS:-24-182:Extrazellulaer 183-205:Transmembran 206-343:Linker 344-435:Todesdomaene}"
SEQ_NEAR_CUTOFF="${SEQ_NEAR_CUTOFF:-8}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#run pipeline
LD_LIBRARY_PATH="${CONDA_PREFIX:-}/lib:${LD_LIBRARY_PATH:-}" python3 \
    "$SCRIPT_DIR/wt_vs_mutant.py" \
    --wt "$WT_PDB" \
    --mut "$MUT_PDB" \
    --chain "$CHAIN" \
    --position "$POSITION" \
    --radii $RADII \
    --domains $DOMAINS \
    --seq-near-cutoff "$SEQ_NEAR_CUTOFF" \
    --out "wt_vs_mut_${OUT_NAME}.csv"

echo ""
echo "Saved."