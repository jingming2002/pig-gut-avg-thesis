#!/bin/bash
# build 3 regulatory trees: lexA, dicA, sinR
# dicA is large, size-controlled; optimized params avoid 17h stall on regulatory
# Usage: nohup bash build_reg_trees.sh > build_reg_trees.log 2>&1 &

cd /lustre/BIF/nobackup/wang500/thesis

VIRUS_DIR="phagcn3_high_apg/gene_families_v2"
BLAST_DIR="phagcn3_high_apg/bact_refs_v2"
CLEAN_BACT="genomad_merged/clean_bacterial_proteins.faa"
WORK="phagcn3_high_apg/trees_v2"

# gene: name | CD-HIT threshold | max bacterial refs (0=unlimited)
declare -A GENES
GENES["lexA_SOS_lysogeny_repressor_K01356"]="0.9 0"
GENES["dicA_regulatory_K22300"]="0.8 400"
GENES["sinR_regulatory_K19449"]="0.9 0"

for gene in lexA_SOS_lysogeny_repressor_K01356 dicA_regulatory_K22300 sinR_regulatory_K19449; do
    params=(${GENES[$gene]})
    CDHIT_C=${params[0]}
    MAX_BACT=${params[1]}
    
    echo ""
    echo "╔══════════════════════════════════════════════════════════╗"
    echo "║ $gene | CD-HIT=$CDHIT_C maxBact=$MAX_BACT | $(date)"
    echo "╚══════════════════════════════════════════════════════════╝"
    
    G="$WORK/$gene"
    mkdir -p "$G"
    VIRUS_FAA="$VIRUS_DIR/${gene}.faa"
    BLAST="$BLAST_DIR/${gene}.tsv"
    VN=$(grep -c "^>" "$VIRUS_FAA")
    echo "  viral sequences: $VN"
    
    # 1. Extract bacterial references (optional downsampling)
    echo "  [1/6] Extract bacterial references..."
    cut -f2 "$BLAST" | sort -u > "$G/bact_ids_all.txt"
    BN_all=$(wc -l < "$G/bact_ids_all.txt")
    if [ "$MAX_BACT" -gt 0 ] && [ "$BN_all" -gt "$MAX_BACT" ]; then
        # randomly downsample to MAX_BACT
        shuf "$G/bact_ids_all.txt" | head -n "$MAX_BACT" > "$G/bact_ids.txt"
        echo "        bacterial refs: $BN_all → downsampled to $MAX_BACT"
    else
        cp "$G/bact_ids_all.txt" "$G/bact_ids.txt"
        echo "        bacterial refs: $BN_all (no downsampling)"
    fi
    seqtk subseq "$CLEAN_BACT" "$G/bact_ids.txt" > "$G/bact_refs.faa"
    
    # 2. Merge + rename
    echo "  [2/6] Merge and label..."
    awk '/^>/{id=substr($1,2); gsub(/[|@].*/,"",id); cnt[id]++; print ">VIRUS_" id "_" cnt[id]; next} {print}' "$VIRUS_FAA" > "$G/combined.faa"
    awk '/^>/{id=substr($1,2); split(id,a,"__"); cnt[a[1]]++; print ">BACT_" a[1] "_" cnt[a[1]]; next} {print}' "$G/bact_refs.faa" >> "$G/combined.faa"
    
    # 3. CD-HIT (dicA uses stricter 0.8)
    echo "  [3/6] CD-HIT (c=$CDHIT_C)..."
    cd-hit -i "$G/combined.faa" -o "$G/nr.faa" -c "$CDHIT_C" -n 5 -M 16000 -d 0 -T 96 > /dev/null 2>&1
    NR=$(grep -c "^>" "$G/nr.faa")
    NV=$(grep -c "^>VIRUS_" "$G/nr.faa")
    NB=$(grep -c "^>BACT_" "$G/nr.faa")
    echo "        after dereplication: $NR (virus$NV + bacteria$NB)"
    
    # 4. MAFFT
    echo "  [4/6] MAFFT..."
    mafft --auto --thread 96 "$G/nr.faa" > "$G/align.fasta" 2>/dev/null
    
    # 5. trimAl
    echo "  [5/6] trimAl..."
    trimal -in "$G/align.fasta" -out "$G/trim.fasta" -gt 0.8 -cons 60 2>/dev/null
    
    # 6. IQ-TREE (optimized: LG+G, UFBoot, avoids slow regulatory runs)
    echo "  [6/6] IQ-TREE (LG+G, UFBoot 1000)..."
    iqtree -s "$G/trim.fasta" -m LG+G -bb 1000 -nt AUTO -ntmax 96 -pre "$G/tree" -redo > /dev/null 2>&1
    
    if [ -f "$G/tree.treefile" ]; then
        echo "        ✓ Done: $NR leaves (virus$NV+bacteria$NB)"
    else
        echo "        ✗ failed"
    fi
    echo "  end: $(date)"
done

echo ""
echo "=== 3 regulatory treesDone ==="
for gene in lexA_SOS_lysogeny_repressor_K01356 dicA_regulatory_K22300 sinR_regulatory_K19449; do
    t="$WORK/$gene/tree.treefile"
    [ -f "$t" ] && echo "  ✓ ${gene%%_*}" || echo "  ✗ ${gene%%_*}"
done
