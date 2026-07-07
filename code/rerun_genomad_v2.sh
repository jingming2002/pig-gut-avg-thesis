#!/bin/bash
# Option B: re-run geNomad (179 dereplicated genomes, incl. new Bacteroidetes)
# merge all genomes -> single geNomad run -> extract clean proteins
# header format as before: >Genus__accession@@contig
cd /lustre/BIF/nobackup/wang500/thesis

GENOMEDIR="ref_genomes"
GENOMAD_DB="PhaGCN3/genomad_db"
MERGED_DIR="genomad_merged_v2"      # new output dir (does not overwrite old)
mkdir -p "$MERGED_DIR"

MERGED_FNA="$MERGED_DIR/all_refs.fna"

echo "###################################################################"
echo "# optionB: re-run geNomad (deduplicateaftergenomesset)"
echo "# start: $(date)"
echo "###################################################################"

# ========== 1. Merge all genomes; prefix contig headers with genome name ==========
echo ""
echo "[1/3] merge $(ls $GENOMEDIR/*.fna | wc -l) genomes..."
> "$MERGED_FNA"
for fna in "$GENOMEDIR"/*.fna; do
    name=$(basename "$fna" .fna)   # Genus__accession
    # contig header: >origContigID -> >Genus__accession@@origContigID
    sed "s/^>/>${name}@@/" "$fna" >> "$MERGED_FNA"
done

n_contigs=$(grep -c '^>' "$MERGED_FNA")
echo "      merged contig count: $n_contigs"
echo "      file size: $(du -sh $MERGED_FNA | cut -f1)"

# ========== 2. Run geNomad (once, all) ==========
echo ""
echo "[2/3] run geNomad (threads=96)..."
echo "      start: $(date)"

genomad end-to-end \
    --cleanup \
    --threads 96 \
    "$MERGED_FNA" \
    "$MERGED_DIR/output" \
    "$GENOMAD_DB" \
    > "$MERGED_DIR/genomad.log" 2>&1

echo "      Done: $(date)"

# ========== 3. Check output ==========
echo ""
echo "[3/3] Check output..."
PROVIRUS_TSV="$MERGED_DIR/output/all_refs_find_proviruses/all_refs_provirus.tsv"
ALL_PROTEINS="$MERGED_DIR/output/all_refs_annotate/all_refs_proteins.faa"

if [ -f "$PROVIRUS_TSV" ]; then
    n_provirus=$(tail -n +2 "$PROVIRUS_TSV" | wc -l)
    echo "      provirus count: $n_provirus"
else
    echo "      ⚠ provirus.tsv not found! check $MERGED_DIR/genomad.log"
fi

if [ -f "$ALL_PROTEINS" ]; then
    n_prot=$(grep -c '^>' "$ALL_PROTEINS")
    echo "      total protein count: $n_prot"
fi

echo ""
echo "###################################################################"
echo "# geNomad Done! $(date)"
echo "# next step: run extract_clean_bacterial_v2.py extractcleanprotein"
echo "###################################################################"
