#!/bin/bash
# count only, no download - how many per genus under complete+RefSeq
cd /lustre/BIF/nobackup/wang500/thesis

GENERA=(
    "Parabacteroides" "Alistipes" "Phocaeicola" "Alloprevotella"
    "Segatella" "Prevotellamassilia" "Xylanibacter" "Sodaliphilus"
)

echo "=== complete + RefSeq criterion(consistent) ==="
for genus in "${GENERA[@]}"; do
    n=$(datasets summary genome taxon "$genus" \
        --assembly-level complete \
        --assembly-source RefSeq \
        --exclude-atypical \
        --as-json-lines 2>/dev/null \
        | grep -o '"accession":"GC[AF]_[0-9]*\.[0-9]*"' \
        | sort -u | wc -l)
    printf "  %-22s %d\n" "$genus" "$n"
done

echo ""
echo "=== relax: complete+chromosome, incl.GenBank(on) ==="
for genus in "${GENERA[@]}"; do
    n=$(datasets summary genome taxon "$genus" \
        --assembly-level complete,chromosome \
        --exclude-atypical \
        --as-json-lines 2>/dev/null \
        | grep -o '"accession":"GC[AF]_[0-9]*\.[0-9]*"' \
        | sort -u | wc -l)
    printf "  %-22s %d\n" "$genus" "$n"
done
