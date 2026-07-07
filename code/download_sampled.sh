#!/bin/bash
cd /lustre/BIF/nobackup/wang500/thesis/ref_genomes

rm -f *.zip *_all.txt *_sample.txt test_acc.txt test.zip

GENERA=(
    "Lactobacillus" "Clostridium" "Streptococcus" "Roseburia"
    "Faecalibacterium" "Blautia" "Ruminococcus"
    "Prevotella" "Bacteroides"
    "Escherichia"
    "Bifidobacterium"
)

N=10

echo "### genusrandomly sample $N genomes | $(date)"

for genus in "${GENERA[@]}"; do
    echo ""
    echo ">>> $genus"
    
    # extract accessions; keep only GC-prefixed (GCF/GCA) with version dot
    datasets summary genome taxon "$genus" \
        --assembly-level complete \
        --assembly-source RefSeq \
        --exclude-atypical \
        --as-json-lines 2>/dev/null \
        | grep -o '"accession":"GC[AF]_[0-9]*\.[0-9]*"' \
        | cut -d'"' -f4 \
        | sort -u > "${genus}_all.txt"
    
    total=$(wc -l < "${genus}_all.txt")
    echo "    available: $total "
    
    if [ "$total" -eq 0 ]; then
        echo "    ⚠ skip"; rm -f "${genus}_all.txt"; continue
    fi
    
    shuf "${genus}_all.txt" | head -n $N > "${genus}_sample.txt"
    
    # download, retry up to 3 times on failure
    for attempt in 1 2 3; do
        datasets download genome accession \
            --inputfile "${genus}_sample.txt" \
            --include genome \
            --filename "${genus}.zip" \
            --no-progressbar 2>&1 | tail -1
        
        # verify zip integrity
        if unzip -t "${genus}.zip" >/dev/null 2>&1; then
            break
        else
            echo "    ⚠  $attempt th download corrupt, retrying..."
            rm -f "${genus}.zip"
            sleep 5
        fi
    done
    
    if [ -f "${genus}.zip" ]; then
        unzip -o -q "${genus}.zip" -d "${genus}_tmp"
        find "${genus}_tmp" -name "*.fna" | while read fna; do
            acc=$(basename "$(dirname "$fna")")
            cp "$fna" "${genus}__${acc}.fna"
        done
        rm -rf "${genus}_tmp" "${genus}.zip" "${genus}_all.txt" "${genus}_sample.txt"
        count=$(ls ${genus}__*.fna 2>/dev/null | wc -l)
        echo "    ✓ download $count "
    else
        echo "    ✗ all threefailed"
        rm -f "${genus}_all.txt" "${genus}_sample.txt"
    fi
    
    sleep 3  # rate limit,avoid NCBI ban
done

echo ""
echo "### Done! $(date)"
echo "total count: $(ls *.fna 2>/dev/null | wc -l)"
echo ""
for genus in "${GENERA[@]}"; do
    printf "  %-20s %d\n" "$genus" "$(ls ${genus}__*.fna 2>/dev/null | wc -l)"
done
echo ""
echo "total size: $(du -sh . 2>/dev/null | cut -f1)"
