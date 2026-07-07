#!/bin/bash
cd /lustre/BIF/nobackup/wang500/thesis
mkdir -p bact_refs_v3

for gene in cysH_sulfur_metabolism_K00390 phoH_phosphate_metabolism_K06217 \
            folD_folate_C1_metabolism_K01491 iscS_cysteine_desulfurase_K04487 \
            mazF_toxin_antitoxin_K07171 spoVT_sporulation_regulation_K04769 \
            lexA_SOS_lysogeny_repressor_K01356 dicA_regulatory_K22300 \
            sinR_regulatory_K19449; do
    echo ">>> $gene | $(date)"
    ./diamond blastp \
        --query phagcn3_high_apg/gene_families_v2/${gene}.faa \
        --db genomad_merged_v2/clean_bact_v2 \
        --out bact_refs_v3/${gene}.tsv \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --threads 96 --max-target-seqs 50 --evalue 1e-5 2>/dev/null
    echo "    hits: $(wc -l < bact_refs_v3/${gene}.tsv)"
    echo "    top donor genera:"
    cut -f2 bact_refs_v3/${gene}.tsv | sed 's/__.*//' | sort | uniq -c | sort -rn | head -8 | sed 's/^/      /'
    echo ""
done
echo "### allDone! $(date)"
