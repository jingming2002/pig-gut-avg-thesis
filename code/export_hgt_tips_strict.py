#!/usr/bin/env python3
"""
truly nested HGTvirusleaves(strict criterion):
virusleavesnearest sister cladepredominantlybacteria(>50%) + dual criterion
- virusclustered together(sister isvirus) -> not counted asHGT(ancient cohesion)
- virusembedded amongbacteriain(sister isbacteria) -> truly nested HGT
"""
from ete3 import Tree
import os

ALRT_MIN = 80
UFBOOT_MIN = 95
SISTER_BACT_FRAC = 0.5   # sister cladebacteriafraction threshold

def ltype(n): return "V" if n.startswith("VIRUS_") else "B"

def parse_support(node):
    name = node.name
    if not name or '/' not in name: return None, None
    try:
        parts = name.split('/')
        return float(parts[0]), float(parts[1])
    except: return None, None

def is_strong(node):
    alrt, ufboot = parse_support(node)
    if alrt is None: return False
    return alrt >= ALRT_MIN and ufboot >= UFBOOT_MIN

genes = [
    "cysH_sulfur_metabolism_K00390","phoH_phosphate_metabolism_K06217",
    "folD_folate_C1_metabolism_K01491","iscS_cysteine_desulfurase_K04487",
    "mazF_toxin_antitoxin_K07171","spoVT_sporulation_regulation_K04769",
    "lexA_SOS_lysogeny_repressor_K01356","dicA_regulatory_K22300",
    "sinR_regulatory_K19449",
]

outdir = "phagcn3_high_apg/hgt_tips_v3_strict"
os.makedirs(outdir, exist_ok=True)

print(f"strict criterion: virusnearest sister cladebacteriafraction>{SISTER_BACT_FRAC} + SH-aLRT>={ALRT_MIN} & UFBoot>={UFBOOT_MIN}")
print("="*70)

for gene in genes:
    path = f"phagcn3_high_apg/trees_v3/{gene}/tree.treefile"
    try:
        t = Tree(path, format=1)
        t.set_outgroup(t.get_midpoint_outgroup())
    except Exception as e:
        print(f"{gene}: failed {e}")
        continue
    
    v_leaves = [l for l in t.get_leaves() if ltype(l.name)=="V"]
    hgt_tips = []
    for vl in v_leaves:
        parent = vl.up
        if parent is None: continue
        # sister clade = all leaves under parent except the virus leaf itself
        sisters = [l for l in parent.get_leaves() if l != vl]
        if not sisters: continue
        # fraction of bacteria among sisters
        n_bact = sum(1 for l in sisters if ltype(l.name)=="B")
        bact_frac = n_bact / len(sisters)
        # Strict criterion: sister predominantly bacterial + node passes dual support
        if bact_frac > SISTER_BACT_FRAC and is_strong(parent):
            hgt_tips.append(vl.name)
    
    with open(f"{outdir}/{gene}.txt", "w") as f:
        for tip in hgt_tips:
            f.write(tip + "\n")
    gname = gene.split("_")[0]
    print(f"{gname:6s}: {len(hgt_tips):3d} truly nested HGT / {len(v_leaves):3d} virus")

print("="*70)
print(f"output to {outdir}/")
