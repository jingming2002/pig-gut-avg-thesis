#!/usr/bin/env python3
"""
cleanthree classesAVGhostdistributioncount
- eachviruscount once(hostphylum/family,deduplicatevirus)
- Chi-square test: three AVG classes x host phylum independence
"""
import csv
from collections import defaultdict
from scipy.stats import chi2_contingency

TABLE = "pig_vOTUs_medium_quality_taxonomy_host.txt"
CSV = {
    'AMG':  "phagcn3_high_amg/final_amg_prediction.csv",
    'APG':  "phagcn3_high_apg/final_apg_prediction.csv",
    'AReG': "phagcn3_high_areg/final_areg_prediction.csv",
}

# ---- 1. Read host table: seq_name -> (phylum, family) ----
# columns (0-based): 0=seq_name, 26=p(phylum), 29=f(family)
host_phylum = {}   # seq_name -> phylum
host_family = {}   # seq_name -> family
with open(TABLE) as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    for row in reader:
        if len(row) < 31:
            continue
        seqname = row[0]
        p = row[26] if len(row)>26 else ""
        fam = row[29] if len(row)>29 else ""
        if p.startswith("p__"):
            host_phylum[seqname] = p
        if fam.startswith("f__"):
            host_family[seqname] = fam

print(f"host: {len(host_phylum)} havephylumprediction, {len(host_family)} havefamilyprediction")

# ---- 2. Per-class AVG viruses -> host phylum (each virus counted once) ----
def phylum_group(p):
    # merge Firmicutes subgroups
    if p.startswith("p__Firmicutes"): return "Firmicutes"
    if p.startswith("p__Bacteroidetes"): return "Bacteroidetes"
    if p.startswith("p__Pseudomonadota"): return "Proteobacteria"
    return "Other"

class_phylum_counts = {}  # cls -> {phylum_group: count}
for cls, csvpath in CSV.items():
    votus = set()
    with open(csvpath) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if row: votus.add(row[0])
    # take each virus's host phylum (count once)
    counts = defaultdict(int)
    n_with_host = 0
    for v in votus:
        if v in host_phylum:
            counts[phylum_group(host_phylum[v])] += 1
            n_with_host += 1
    class_phylum_counts[cls] = counts
    print(f"{cls}: {len(votus)} virus, {n_with_host} hostphylumprediction")

# ---- 3. Build contingency table (phylum level) ----
phyla = ["Firmicutes","Bacteroidetes","Proteobacteria","Other"]
print("\n=== contingency table (viruscount, virusonce) ===")
print(f"{'':6s}", "  ".join(f"{p[:12]:>13s}" for p in phyla))
table = []
for cls in ['AMG','APG','AReG']:
    row = [class_phylum_counts[cls].get(p,0) for p in phyla]
    table.append(row)
    print(f"{cls:6s}", "  ".join(f"{v:>13d}" for v in row))

# ---- 4. Chi-square test ----
import numpy as np
table = np.array(table)
chi2, p, dof, expected = chi2_contingency(table)
print(f"\n=== Chi-square test (three classesAVG x hostphylum) ===")
print(f"Chi2 = {chi2:.3f}, dof = {dof}, p = {p:.4e}")
print(f": {expected.min():.1f} (>5 chi-square is reliable)")
if p < 0.05:
    print(">>> p < 0.05: three classesAVGhostphylumdistributioncount")
else:
    print(">>> p >= 0.05: countthree classeshostphylumdistribution")

# ---- 5. Per-class proportions (for interpretation) ----
print("\n=== hostphylumfraction ===")
for i,cls in enumerate(['AMG','APG','AReG']):
    tot = table[i].sum()
    print(f"{cls:6s}", "  ".join(f"{p[:6]}:{100*table[i][j]/tot:.0f}%" for j,p in enumerate(phyla)))
