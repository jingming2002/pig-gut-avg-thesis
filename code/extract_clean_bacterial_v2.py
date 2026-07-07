#!/usr/bin/env python3
"""
based on provirus coordinates,fromallbacteriaexclude proteins falling within provirus genes in region,
obtain cleanbacteriaprotein。
method basis: Baeuerle et al. 2026 - use geNomad identify and exclude prophage contamination
"""
import re
from collections import defaultdict

PROVIRUS_TSV = "genomad_merged_v2/output/all_refs_find_proviruses/all_refs_provirus.tsv"
ALL_PROTEINS = "genomad_merged_v2/output/all_refs_annotate/all_refs_proteins.faa"
CLEAN_OUT = "genomad_merged_v2/clean_bacterial_proteins.faa"

# ============================================================
# 1. Read provirus coordinates: each contig -> [(start, end), ...]
# ============================================================
print("[1/3] read provirus coordinates...")
provirus_regions = defaultdict(list)
with open(PROVIRUS_TSV) as f:
    next(f)  # skip header
    for line in f:
        cols = line.rstrip('\n').split('\t')
        source_seq = cols[1]   # source contig
        start = int(cols[2])
        end = int(cols[3])
        provirus_regions[source_seq].append((start, end))

n_contigs_with_provirus = len(provirus_regions)
n_provirus = sum(len(v) for v in provirus_regions.values())
print(f"      {n_provirus}  provirus, distribution {n_contigs_with_provirus}  contig on")

# ============================================================
# 2. Iterate all proteins, parse coordinates, exclude those within provirus regions
# ============================================================
print("[2/3] exclude provirus genes in region...")

# all-protein header format:
# >Agathobacter__GCA_000020605.1@@CP001107.1_1 # 1 # 1362 # 1 # ID=...
# contig = between @@ and last _; gene coordinates at # start # end #

def parse_header(header):
    """Return (contig, gene_start, gene_end)"""
    # strip > from header
    parts = header[1:].split(' # ')
    gene_id = parts[0]          # Agathobacter__...CP001107.1_1
    start = int(parts[1])
    end = int(parts[2])
    # contig = gene_id without trailing _number
    contig = gene_id.rsplit('_', 1)[0]
    return gene_id, contig, start, end

def in_provirus(contig, gstart, gend):
    """Whether gene coordinates fall within any provirus region of the contig"""
    if contig not in provirus_regions:
        return False
    # gene midpoint within region counts (overlap also possible)
    mid = (gstart + gend) / 2
    for (ps, pe) in provirus_regions[contig]:
        if ps <= mid <= pe:
            return True
    return False

total = 0
excluded = 0
kept = 0

with open(ALL_PROTEINS) as inf, open(CLEAN_OUT, 'w') as outf:
    header = None
    seq_lines = []
    
    def flush():
        global excluded, kept
        if header is None:
            return
        gene_id, contig, gs, ge = parse_header(header)
        if in_provirus(contig, gs, ge):
            excluded += 1
        else:
            outf.write(header + '\n')
            outf.write(''.join(seq_lines))
            kept += 1
    
    for line in inf:
        if line.startswith('>'):
            flush()
            header = line.rstrip('\n')
            seq_lines = []
            total += 1
        else:
            seq_lines.append(line)
    flush()  # after

print(f"      allprotein: {total}")
print(f"      exclude(provirus): {excluded}")
print(f"      keep(cleanbacteria): {kept}")

# ============================================================
# 3. Done
# ============================================================
print(f"[3/3] Done! cleanbacteriaprotein -> {CLEAN_OUT}")
print(f"\nnext step: use 8  AVG sequence blastp search this clean database")
