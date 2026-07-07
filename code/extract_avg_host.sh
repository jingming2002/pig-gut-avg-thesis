#!/bin/bash
# extract rows of AVG-carrying viruses from CheckV+iPHoP table (completeness + host prediction)
cd /lustre/BIF/nobackup/wang500/thesis

TABLE="pig_vOTUs_medium_quality_taxonomy_host.txt"
GENEDIR="phagcn3_high_apg/gene_families_v2"
OUT="avg_virus_host"
mkdir -p "$OUT"

genes=(
    "cysH_sulfur_metabolism_K00390"
    "phoH_phosphate_metabolism_K06217"
    "folD_folate_C1_metabolism_K01491"
    "iscS_cysteine_desulfurase_K04487"
    "mazF_toxin_antitoxin_K07171"
    "spoVT_sporulation_regulation_K04769"
    "lexA_SOS_lysogeny_repressor_K01356"
    "dicA_regulatory_K22300"
    "sinR_regulatory_K19449"
)

# header
head -1 "$TABLE" > "$OUT/header.txt"

# collect contig IDs of all AVG viruses (part before @)
> /tmp/all_avg_ids.txt
echo "=== per geneAVGvirusmatch ==="
for gene in "${genes[@]}"; do
    grep "^>" "$GENEDIR/${gene}.faa" | sed 's/^>//; s/@.*//' | sort -u > /tmp/${gene}_ids.txt
    n_virus=$(wc -l < /tmp/${gene}_ids.txt)
    cat /tmp/${gene}_ids.txt >> /tmp/all_avg_ids.txt
    
    # rows of this gene's viruses in table (matched by part before @)
    awk 'NR==FNR{ids[$1]=1; next} FNR==1{next} {key=$1; sub(/@.*/,"",key); if(key in ids) print}' \
        /tmp/${gene}_ids.txt "$TABLE" > "$OUT/${gene}_host.tsv"
    n_matched=$(wc -l < "$OUT/${gene}_host.tsv")
    # with high-confidence host prediction (Confidence non-empty; roughly via Host column)
    printf "  %-40s virus %3d | inmatch %3d\n" "$gene" "$n_virus" "$n_matched"
done

# merge table rows of all AVG viruses (deduplicated)
sort -u /tmp/all_avg_ids.txt > /tmp/all_avg_ids_uniq.txt
echo ""
echo "allAVGvirus contig(deduplicate): $(wc -l < /tmp/all_avg_ids_uniq.txt)"

cat "$OUT/header.txt" > "$OUT/all_avg_virus_host.tsv"
awk 'NR==FNR{ids[$1]=1; next} FNR==1{next} {key=$1; sub(/@.*/,"",key); if(key in ids) print}' \
    /tmp/all_avg_ids_uniq.txt "$TABLE" >> "$OUT/all_avg_virus_host.tsv"
n_total=$(tail -n +2 "$OUT/all_avg_virus_host.tsv" | wc -l)
echo "matchAVGvirus: $n_total"

# with host prediction (Host column 24 non-empty)
echo ""
echo "=== host prediction summary ==="
n_host=$(tail -n +2 "$OUT/all_avg_virus_host.tsv" | awk -F'\t' '$24!="" {c++} END{print c+0}')
echo "with host predictionAVGvirus: $n_host"

echo ""
echo "=== hostphylumdistribution(with prediction) ==="
tail -n +2 "$OUT/all_avg_virus_host.tsv" | awk -F'\t' '$26!=""{print $26}' | sort | uniq -c | sort -rn | head -15

echo ""
echo "output dir: $OUT/"
echo "  all_avg_virus_host.tsv (all)"
echo "  <gene>_host.tsv (per gene)"
