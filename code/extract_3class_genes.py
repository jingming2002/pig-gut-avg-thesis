#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
import os

BASEDIR = "checkamg_output/results"
OUTDIR = "phagcn3_high_apg/gene_families"
os.makedirs(OUTDIR, exist_ok=True)

# select most representative per class
TARGETS = {
    'K00390': ('cysH_sulfate_metabolism', 'metabolic'),
    'K06217': ('phoH_phosphate_metabolism', 'metabolic'),
    'K07171': ('mazF_toxin_antitoxin', 'physiological'),
    'K07729': ('transcriptional_regulator', 'regulatory'),
    'K03111': ('ssb_DNA_binding', 'regulatory')
}

print("Reading final_results.tsv...")
df = pd.read_csv(f"{BASEDIR}/final_results.tsv", sep="\t", low_memory=False)

print("Loading high-confidence AVG sequences...")
fasta_files = [
    f"{BASEDIR}/faa_metabolic/AMGs_high_confidence.faa",
    f"{BASEDIR}/faa_physiology/APGs_high_confidence.faa",
    f"{BASEDIR}/faa_regulatory/AReGs_high_confidence.faa"
]

all_seqs = {}
for fasta in fasta_files:
    if os.path.exists(fasta):
        for rec in SeqIO.parse(fasta, 'fasta'):
            all_seqs[rec.id] = rec
        print(f"  ✓ {os.path.basename(fasta)}")

print(f"Total: {len(all_seqs)} sequences\n")

print("Extracting 5 representative genes from 3 AVG classes:")
print("="*90)

summary = {}
for ko, (name, cls) in TARGETS.items():
    proteins = df[(df['KEGG KO'] == ko) & (df['Protein Viral Origin Confidence'] == 'high')]['Protein'].dropna().unique()
    
    found_count = 0
    out_file = f"{OUTDIR}/{name}_{ko}.faa"
    
    with open(out_file, 'w') as out:
        for prot_id in proteins:
            if prot_id in all_seqs:
                SeqIO.write(all_seqs[prot_id], out, 'fasta')
                found_count += 1
    
    summary[ko] = found_count
    status = "✓" if found_count > 0 else "✗"
    print(f"{status} {cls:15s} | {ko} ({name:35s}): {found_count:5d} sequences")

print("="*90)
print(f"Total extracted: {sum(summary.values())} sequences\n")
os.system(f"ls -lh {OUTDIR}/")

