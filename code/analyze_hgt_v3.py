#!/usr/bin/env python3
"""
Qualitative HGT analysis v3 (dual support: SH-aLRT>=80 and UFBoot>=95)
- strongly-supportedusevirusmonophyly/HGT
- new reference with balanced Bacteroidetes
"""
from ete3 import Tree
import statistics

# dual support-value thresholds
ALRT_MIN = 80
UFBOOT_MIN = 95

PHYLUM = {
    'Clostridium':'Firmicutes','Ruminococcus':'Firmicutes','Blautia':'Firmicutes',
    'Faecalibacterium':'Firmicutes','Streptococcus':'Firmicutes','Lactobacillus':'Firmicutes',
    'Roseburia':'Firmicutes','Oscillibacter':'Firmicutes','Agathobacter':'Firmicutes',
    'Megasphaera':'Firmicutes','Phascolarctobacterium':'Firmicutes','Eubacterium':'Firmicutes',
    'Bacteroides':'Bacteroidetes','Prevotella':'Bacteroidetes','Parabacteroides':'Bacteroidetes',
    'Alistipes':'Bacteroidetes','Phocaeicola':'Bacteroidetes','Segatella':'Bacteroidetes',
    'Treponema':'Spirochaetes','Bifidobacterium':'Actinobacteria','Escherichia':'Proteobacteria',
}
def phylum_of(g): return PHYLUM.get(g,'Other')
def ltype(n): return "V" if n.startswith("VIRUS_") else "B"
def genus_of(n):
    parts = n.split("_")
    return parts[1] if len(parts)>1 else "?"

def parse_support(node):
    """Parse node SH-aLRT/UFBoot; return (alrt, ufboot) or (None, None)"""
    name = node.name
    if not name or '/' not in name:
        return None, None
    try:
        parts = name.split('/')
        alrt = float(parts[0])
        ufboot = float(parts[1])
        return alrt, ufboot
    except:
        return None, None

def is_strong(node):
    """Dual criterion: SH-aLRT>=80 and UFBoot>=95"""
    alrt, ufboot = parse_support(node)
    if alrt is None: return False
    return alrt >= ALRT_MIN and ufboot >= UFBOOT_MIN

genes = [
    ("cysH","sulfur_metabolism_K00390","Metabolic"),
    ("phoH","phosphate_metabolism_K06217","Metabolic"),
    ("folD","folate_C1_metabolism_K01491","Metabolic"),
    ("iscS","cysteine_desulfurase_K04487","Metabolic"),
    ("mazF","toxin_antitoxin_K07171","Physiological"),
    ("spoVT","sporulation_regulation_K04769","Physiological"),
    ("lexA","SOS_lysogeny_repressor_K01356","Regulatory"),
    ("dicA","regulatory_K22300","Regulatory"),
    ("sinR","regulatory_K19449","Regulatory"),
]

def analyze(path):
    t = Tree(path, format=1)
    try: t.set_outgroup(t.get_midpoint_outgroup())
    except: pass
    leaves = t.get_leaves()
    v_leaves = [l for l in leaves if ltype(l.name)=="V"]
    nv = len(v_leaves)
    nb = sum(1 for l in leaves if ltype(l.name)=="B")
    
    # largest pure-viral cluster (strongly-supported clades only)
    pure_v_strong = []
    for node in t.traverse():
        if node.is_leaf(): continue
        lv = node.get_leaves()
        if all(ltype(l.name)=="V" for l in lv):
            # only strongly-supported pure-viral clusters count
            if is_strong(node):
                pure_v_strong.append(len(lv))
    max_v = max(pure_v_strong) if pure_v_strong else 1
    frac = max_v/nv if nv>0 else 0
    
    # number of viral clusters (strongly-supported independent pure-viral clusters)
    v_clusters = 0; visited=set()
    for node in sorted([n for n in t.traverse() if not n.is_leaf()],
                       key=lambda x:-len(x.get_leaves())):
        lv=node.get_leaves()
        if all(ltype(l.name)=="V" for l in lv) and is_strong(node):
            ids=frozenset(l.name for l in lv)
            if not any(l in visited for l in ids):
                v_clusters+=1; visited.update(ids)
    for vl in v_leaves:
        if vl.name not in visited: v_clusters+=1; visited.add(vl.name)
    
    # HGT event: virus nested in bacterial clade, node strongly supported
    # (virus sister is bacterial, and parent node strongly supported)
    hgt_events = 0
    donor_genera = {}
    donor_phyla = {}
    for vl in v_leaves:
        node = vl
        while node.up:
            parent = node.up
            sibs = [l for l in parent.get_leaves() if l != vl]
            has_bact = any(ltype(l.name)=="B" for l in sibs)
            if has_bact:
                # find nearest bacteria-containing ancestor, check node support
                if is_strong(parent):
                    hgt_events += 1
                    # take only the single nearest bacterium as donor (avoids large-clade inflation)
                    bact_sibs = [l for l in parent.get_leaves() if ltype(l.name)=="B"]
                    nearest = min(bact_sibs, key=lambda b: vl.get_distance(b))
                    g = genus_of(nearest.name)
                    donor_genera[g] = donor_genera.get(g,0)+1
                    ph = phylum_of(g)
                    donor_phyla[ph] = donor_phyla.get(ph,0)+1
                break
            node = parent
    
    return {'nv':nv,'nb':nb,'frac':frac,'v_clusters':v_clusters,
            'hgt_events':hgt_events,'donor_genera':donor_genera,'donor_phyla':donor_phyla}

print("="*100)
print(f"HGTqualitative analysis v3 (dual criterion: SH-aLRT>={ALRT_MIN} and UFBoot>={UFBOOT_MIN})")
print("with balanced Bacteroidetes reference")
print("="*100)
print(f"{'Gene':6s} {'Class':12s} {'V':>4s} {'B':>4s} {'MaxVclust%':>10s} {'#Vclust':>7s} {'HGT(strong)':>11s}  {'phylumdonor':30s}")
print("-"*100)

results=[]
for short, suffix, cls in genes:
    path = f"phagcn3_high_apg/trees_v3/{short}_{suffix}/tree.treefile"
    try:
        r = analyze(path)
    except Exception as e:
        print(f"{short:6s} failed: {e}")
        continue
    # phylum donor summary
    ph_sum = ", ".join(f"{p}:{c}" for p,c in sorted(r['donor_phyla'].items(),key=lambda x:-x[1])[:3])
    print(f"{short:6s} {cls:12s} {r['nv']:4d} {r['nb']:4d} {100*r['frac']:9.1f}% {r['v_clusters']:7d} {r['hgt_events']:11d}  {ph_sum}")
    results.append((short,cls,r))

# three-class summary
print("\n"+"="*100)
print("three-class summary")
print("="*100)
for cls in ['Metabolic','Physiological','Regulatory']:
    sub=[r for r in results if r[1]==cls]
    if not sub: continue
    fracs=[r[2]['frac'] for r in sub]
    clusts=[r[2]['v_clusters'] for r in sub]
    hgts=[r[2]['hgt_events'] for r in sub]
    # merge phylum donors for this class
    ph_total={}
    for r in sub:
        for p,c in r[2]['donor_phyla'].items():
            ph_total[p]=ph_total.get(p,0)+c
    print(f"\n{cls} ({len(sub)}gene):")
    print(f"  meanviruscluster: {100*statistics.mean(fracs):.1f}%")
    print(f"  meannumber of viral clusters: {statistics.mean(clusts):.1f}")
    print(f"  total strongly-supported HGT events: {sum(hgts)}")
    print(f"  phylumdonor distribution: " + ", ".join(f"{p}:{c}" for p,c in sorted(ph_total.items(),key=lambda x:-x[1])))

# Bacteroidetes involvement comparison (key)
print("\n"+"="*100)
print("⭐ Bacteroidetes involvement (validate reference bias)")
print("="*100)
for short,cls,r in results:
    bact = r['donor_phyla'].get('Bacteroidetes',0)
    firm = r['donor_phyla'].get('Firmicutes',0)
    total = sum(r['donor_phyla'].values())
    if total>0:
        print(f"  {short:6s} ({cls:12s}): Bacteroidetes {bact:3d} / Firmicutes {firm:3d} "
              f"(Bactfraction {100*bact/total:.0f}%)")
