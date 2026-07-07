#!/bin/bash
cd /lustre/BIF/nobackup/wang500/thesis

OUTDIR="ref_genomes"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# pig gut core genera (11)
GENERA=(
    "Lactobacillus" "Clostridium" "Streptococcus" "Roseburia"
    "Faecalibacterium" "Blautia" "Ruminococcus"
    "Prevotella" "Bacteroides"
    "Escherichia"
    "Bifidobacterium"
)

echo "### downloadpig gut core microbiome genomes | start: $(date)"
echo "### strategy(tier3): genus 20 genomes (complete + chromosome)"

for genus in "${GENERA[@]}"; do
    echo ""
    echo ">>> $genus"
    
    datasets download genome taxon "$genus" \
        --assembly-level complete,chromosome \
        --include genome \
        --limit 20 \
        --filename "${genus}.zip" \
        --no-progressbar 2>&1 | tail -2
    
    if [ -f "${genus}.zip" ]; then
        unzip -o -q "${genus}.zip" -d "${genus}_tmp"
        find "${genus}_tmp" -name "*.fna" | while read fna; do
            acc=$(basename "$(dirname "$fna")")
            cp "$fna" "${genus}__${acc}.fna"
        done
        rm -rf "${genus}_tmp" "${genus}.zip"
        count=$(ls ${genus}__*.fna 2>/dev/null | wc -l)
        echo "    ✓ $count genomes"
    else
        echo "    ✗ downloadfailed"
    fi
done

echo ""
echo "### Done! $(date)"
echo "total genomes: $(ls *.fna 2>/dev/null | wc -l)"
echo ""
echo "genusdistribution:"
for genus in "${GENERA[@]}"; do
    count=$(ls ${genus}__*.fna 2>/dev/null | wc -l)
    printf "  %-20s %d\n" "$genus" "$count"
done
echo ""
echo "total size: $(du -sh . 2>/dev/null | cut -f1)"
