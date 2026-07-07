#!/usr/bin/env python3
"""
from clean AVG Dataextract 8 tree buildinggeneproteinsequence
by KEGG KO match,from three classes FASTA fetch sequences from file
"""
import pandas as pd
import os

# 8 target genes
TARGETS = {
    "K00390": ("cysH",   "sulfur_metabolism",      "metabolic"),
    "K06217": ("phoH",   "phosphate_metabolism",   "metabolic"),
    "K01491": ("folD",   "folate_C1_metabolism",   "metabolic"),
    "K04487": ("iscS",   "cysteine_desulfurase",   "metabolic"),
    "K07171": ("mazF",   "toxin_antitoxin",        "physiological"),
    "K04769": ("spoVT",  "sporulation_regulation", "physiological"),
    "K07729": ("reg",    "transcriptional_reg",    "regulatory"),
    "K01356": ("lexA",   "SOS_lysogeny_repressor", "regulatory"),
}

# three-class FASTA paths
FAA = {
    "metabolic":     "checkamg_output/results/faa_metabolic/AMGs_high_confidence.faa",
    "physiological": "checkamg_output/results/faa_physiology/APGs_high_confidence.faa",
    "regulatory":    "checkamg_output/results/faa_regulatory/AReGs_high_confidence.faa",
}

OUTDIR = "phagcn3_high_apg/gene_families_v2"
os.makedirs(OUTDIR, exist_ok=True)

# read final_results, build Protein -> KEGG KO mapping
print("read final_results.tsv...")
df = pd.read_csv("checkamg_output/results/final_results.tsv", sep="\t", low_memory=False)
high = df[df['Protein Viral Origin Confidence'] == 'high']

# set of Protein IDs for each target KO
ko_proteins = {}
for ko in TARGETS:
    prots = set(high[high['KEGG KO'] == ko]['Protein'].astype(str))
    ko_proteins[ko] = prots
    print(f"  {ko} ({TARGETS[ko][0]}): {len(prots)}  Protein ID")

# helper function to read FASTA
def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                # ID is before first space in header
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = ''.join(buf)
    return seqs

# load three-class FASTA
print("\nload three-class FASTA...")
all_seqs = {}
for cls, path in FAA.items():
    s = read_fasta(path)
    all_seqs[cls] = s
    print(f"  {cls}: {len(s)} sequences")

# extract each target gene
print("\nextract sequences...")
for ko, (gene, desc, cls) in TARGETS.items():
    prots = ko_proteins[ko]
    seqs = all_seqs[cls]
    
    out_path = f"{OUTDIR}/{gene}_{desc}_{ko}.faa"
    written = 0
    with open(out_path, 'w') as out:
        for pid in prots:
            if pid in seqs:
                out.write(f">{pid}\n{seqs[pid]}\n")
                written += 1
    print(f"  {gene} ({ko}): write {written}  → {out_path}")

print(f"\nDone! 8 gene sequences in {OUTDIR}/")
print("\nsequence count per file:")
for ko, (gene, desc, cls) in TARGETS.items():
    p = f"{OUTDIR}/{gene}_{desc}_{ko}.faa"
    if os.path.exists(p):
        n = sum(1 for l in open(p) if l.startswith('>'))
        print(f"  {gene:8s} {n:5d}")
