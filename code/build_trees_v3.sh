#!/bin/bash
# ===================================================================
# 9 AVG gene HGT tree building v3 (new reference with balanced Bacteroidetes)
# Usage: nohup bash build_trees_v3.sh > build_trees_v3.log 2>&1 &
# ===================================================================
cd /lustre/BIF/nobackup/wang500/thesis

VIRUS_DIR="phagcn3_high_apg/gene_families_v2"
BLAST_DIR="bact_refs_v3"                                  # new blast results
CLEAN_BACT="genomad_merged_v2/clean_bacterial_proteins.faa"  # new clean database
WORK="phagcn3_high_apg/trees_v3"                          # new output (does not overwrite v2)
mkdir -p "$WORK"

# 9 genes (dicA marked for downsampling)
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

echo "###################################################################"
echo "# 9 gene HGT tree building v3 (incl.Bacteroidetes) | start: $(date)"
echo "###################################################################"

for gene in "${genes[@]}"; do
    echo ""
    echo "╔══════════════════════════════════════════════════════════╗"
    echo "║ $gene | $(date)"
    echo "╚══════════════════════════════════════════════════════════╝"
    
    G="$WORK/$gene"
    mkdir -p "$G"
    
    VIRUS_FAA="$VIRUS_DIR/${gene}.faa"
    BLAST="$BLAST_DIR/${gene}.tsv"
    
    VN=$(grep -c "^>" "$VIRUS_FAA")
    echo "  viral sequences: $VN"
    
    # ---- 1. Extract bacterial references ----
    echo "  [1/6] Extract bacterial references..."
    cut -f2 "$BLAST" | sort -u > "$G/bact_ids.txt"
    BN_ids=$(wc -l < "$G/bact_ids.txt")
    
    # dicA downsampling: too many bacterial refs, randomly sample 400
    if [ "$gene" = "dicA_regulatory_K22300" ]; then
        echo "        dicA: bacterial refs $BN_ids , downsampled to 400"
        shuf "$G/bact_ids.txt" | head -400 > "$G/bact_ids_sub.txt"
        mv "$G/bact_ids_sub.txt" "$G/bact_ids.txt"
        BN_ids=$(wc -l < "$G/bact_ids.txt")
    fi
    
    seqtk subseq "$CLEAN_BACT" "$G/bact_ids.txt" > "$G/bact_refs.faa"
    BN=$(grep -c "^>" "$G/bact_refs.faa")
    echo "        bacterial refs: $BN"
    
    # ---- 2. Merge and label source ----
    echo "  [2/6] Merge and label source..."
    awk '/^>/{id=substr($1,2); gsub(/[|@].*/,"",id); cnt[id]++;
             print ">VIRUS_" id "_" cnt[id]; next} {print}' "$VIRUS_FAA" > "$G/combined.faa"
    awk '/^>/{id=substr($1,2); split(id,a,"__"); cnt[a[1]]++;
             print ">BACT_" a[1] "_" cnt[a[1]]; next} {print}' "$G/bact_refs.faa" >> "$G/combined.faa"
    TOTAL=$(grep -c "^>" "$G/combined.faa")
    echo "        total: $TOTAL"
    
    # ---- 3. CD-HIT ----
    echo "  [3/6] CD-HIT dereplication..."
    # dicA uses more aggressive c=0.8
    CDHIT_C=0.9
    if [ "$gene" = "dicA_regulatory_K22300" ]; then CDHIT_C=0.8; fi
    cd-hit -i "$G/combined.faa" -o "$G/nr.faa" -c $CDHIT_C -n 5 -M 16000 -d 0 -T 96 > /dev/null 2>&1
    NR=$(grep -c "^>" "$G/nr.faa")
    NV=$(grep -c "^>VIRUS_" "$G/nr.faa")
    NB=$(grep -c "^>BACT_" "$G/nr.faa")
    echo "        after dereplication: $NR (virus $NV + bacteria $NB)"
    
    # ---- 4. MAFFT ----
    echo "  [4/6] MAFFT alignment..."
    mafft --auto --thread 96 "$G/nr.faa" > "$G/align.fasta" 2>/dev/null
    
    # ---- 5. trimAl ----
    echo "  [5/6] trimAl trimming..."
    trimal -in "$G/align.fasta" -out "$G/trim.fasta" -gt 0.8 -cons 60 2>/dev/null
    
    # ---- 6. IQ-TREE (optimized: LG+G, not TEST, avoids stalling) ----
    echo "  [6/6] IQ-TREE..."
    iqtree -s "$G/trim.fasta" -m LG+G -bb 1000 -alrt 1000 -nt AUTO -ntmax 96 -pre "$G/tree" -quiet -redo 2>/dev/null
    
    if [ -f "$G/tree.treefile" ]; then
        echo "        ✓ treeDone: $NR leaves (virus $NV + bacteria $NB)"
    else
        echo "        ✗ tree building failed"
    fi
    echo "  end: $(date)"
done

echo ""
echo "###################################################################"
echo "# allDone! $(date)"
echo "###################################################################"
for gene in "${genes[@]}"; do
    t="$WORK/$gene/tree.treefile"
    if [ -f "$t" ]; then
        leaves=$(grep -o "," "$t" | wc -l)
        echo "  ✓ $gene ($((leaves+1)) leaves)"
    else
        echo "  ✗ $gene (failed)"
    fi
done

