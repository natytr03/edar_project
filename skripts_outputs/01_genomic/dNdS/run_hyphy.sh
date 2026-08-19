#!/usr/bin/env bash


set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"

#for me: activate HyPhy environment and correct directory
conda activate hyphy_env
cd ~/edar/msa/codonmsa


ALIGNMENT="species_list_edar_codon_5.final_aligned.fasta"
TREE="edar_paml_tree_5.final.nwk"

echo "RUN FUBAR"
hyphy FUBAR \
    --alignment "$ALIGNMENT" \
    --tree "$TREE" \
    --output species_list_edar_5.final_fubar.json \
    2>&1 | tee species_list_edar_5.final_fubar.log

echo "RUN FEL"
hyphy FEL \
    --alignment "$ALIGNMENT" \
    --tree "$TREE" \
    --pvalue 0.1 \
    --output species_list_edar_5.final_fel.json \
    2>&1 | tee species_list_edar_5.final_fel.log

echo "RUN MEME"
hyphy MEME \
    --alignment "$ALIGNMENT" \
    --tree "$TREE" \
    --pvalue 0.1 \
    --output species_list_edar_5.final_meme.json \
    2>&1 | tee species_list_edar_5.final_meme.log

echo "RUN aBSREL"
hyphy aBSREL \
    --alignment "$ALIGNMENT" \
    --tree "$TREE" \
    --output species_list_edar_5.final_absrel.json \
    2>&1 | tee species_list_edar_5.final_absrel.log

echo
echo "All HyPhy analyses completet!"
