#!/usr/bin/env python3
"""
Core HGT figures (field-standard):
figureA: three classesAVGviruscohesion comparison(metabolic cohesion/ancient → regulatory dispersal/recentHGT)
Panel B: donor genus heatmap(showFirmicutesdominant)
recompute directly from trees,ensureDataconsistent
"""
from ete3 import Tree
import statistics
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch

mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

PHYLUM = {
    'Clostridium':'Firmicutes','Ruminococcus':'Firmicutes','Blautia':'Firmicutes',
    'Faecalibacterium':'Firmicutes','Streptococcus':'Firmicutes','Lactobacillus':'Firmicutes',
    'Roseburia':'Firmicutes','Oscillibacter':'Firmicutes','Agathobacter':'Firmicutes',
    'Megasphaera':'Firmicutes','Phascolarctobacterium':'Firmicutes','Eubacterium':'Firmicutes',
    'Coprococcus':'Firmicutes','Dorea':'Firmicutes','Lachnoclostridium':'Firmicutes',
    'Anaerostipes':'Firmicutes','Subdoligranulum':'Firmicutes',
    'Bacteroides':'Bacteroidetes','Prevotella':'Bacteroidetes','Parabacteroides':'Bacteroidetes',
    'Alistipes':'Bacteroidetes','Treponema':'Spirochaetes',
    'Bifidobacterium':'Actinobacteria','Collinsella':'Actinobacteria',
    'Escherichia':'Proteobacteria','Desulfovibrio':'Proteobacteria',
}
def phylum_of(g): return PHYLUM.get(g, 'Other')

def ltype(n): return "V" if n.startswith("VIRUS_") else "B"
def genus_of(n):
    parts = n.split("_")
    return parts[1] if len(parts) > 1 else "?"

genes = [
    ("cysH","sulfur_metabolism_K00390","Metabolic","cysH"),
    ("phoH","phosphate_metabolism_K06217","Metabolic","phoH"),
    ("folD","folate_C1_metabolism_K01491","Metabolic","folD"),
    ("iscS","cysteine_desulfurase_K04487","Metabolic","iscS"),
    ("mazF","toxin_antitoxin_K07171","Physiological","mazF"),
    ("spoVT","sporulation_regulation_K04769","Physiological","spoVT"),
    ("lexA","SOS_lysogeny_repressor_K01356","Regulatory","lexA"),
    ("dicA","regulatory_K22300","Regulatory","dicA"),
    ("sinR","regulatory_K19449","Regulatory","sinR"),
]

CLASS_COLORS = {'Metabolic':'#0077BB','Physiological':'#009E73','Regulatory':'#AA3DDD'}

def analyze_tree(path):
    t = Tree(path, format=1)
    try: t.set_outgroup(t.get_midpoint_outgroup())
    except: pass
    leaves = t.get_leaves()
    v_leaves = [l for l in leaves if ltype(l.name)=="V"]
    nv = len(v_leaves)
    nb = sum(1 for l in leaves if ltype(l.name)=="B")
    
    # largest pure-viral cluster
    pure_v = [len(n.get_leaves()) for n in t.traverse() 
              if not n.is_leaf() and all(ltype(l.name)=="V" for l in n.get_leaves())]
    max_v = max(pure_v) if pure_v else 1
    frac = max_v/nv if nv>0 else 0
    
    # number of viral clusters
    v_clusters = 0; visited=set()
    for node in sorted([n for n in t.traverse() if not n.is_leaf()], key=lambda x:-len(x.get_leaves())):
        lv=node.get_leaves()
        if all(ltype(l.name)=="V" for l in lv):
            ids=frozenset(l.name for l in lv)
            if not any(l in visited for l in ids):
                v_clusters+=1; visited.update(ids)
    for vl in v_leaves:
        if vl.name not in visited: v_clusters+=1; visited.add(vl.name)
    
    # donor genus (nearest neighbour)
    neighbor_genera={}
    for vl in v_leaves:
        node=vl
        while node.up:
            node=node.up
            sibs=[l for l in node.get_leaves() if ltype(l.name)=="B"]
            if sibs:
                for b in sibs:
                    g=genus_of(b.name)
                    neighbor_genera[g]=neighbor_genera.get(g,0)+1
                break
    
    return {'nv':nv,'nb':nb,'frac':frac,'v_clusters':v_clusters,'donors':neighbor_genera}

# collectData
data=[]
for short,suffix,cls,label in genes:
    r=analyze_tree(f"phagcn3_high_apg/trees_v2/{short}_{suffix}/tree.treefile")
    data.append({'gene':label,'cls':cls,**r})

# ============ Panel A: cohesion comparison across three classes ============
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# sort genes by class
order = sorted(range(len(data)), key=lambda i:(['Metabolic','Physiological','Regulatory'].index(data[i]['cls']), -data[i]['frac']))
genes_ord = [data[i] for i in order]
labels = [d['gene'] for d in genes_ord]
fracs = [100*d['frac'] for d in genes_ord]
colors = [CLASS_COLORS[d['cls']] for d in genes_ord]

# left: largest viral cluster fraction (higher = older)
bars1 = ax1.bar(range(len(labels)), fracs, color=colors, edgecolor='black', linewidth=0.6)
ax1.set_xticks(range(len(labels)))
ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=11, fontstyle='italic')
ax1.set_ylabel('Largest viral cluster\n(% of viral sequences)', fontsize=12)
ax1.set_title('Viral sequence cohesion', fontsize=13, fontweight='bold')
ax1.axhline(50, ls='--', color='gray', lw=0.8, alpha=0.6)
ax1.grid(axis='y', alpha=0.3)
# class-mean annotation
# right: number of viral clusters (more = dispersed/recent HGT) - log scale
nclusts = [d['v_clusters'] for d in genes_ord]
bars2 = ax2.bar(range(len(labels)), nclusts, color=colors, edgecolor='black', linewidth=0.6)
ax2.set_xticks(range(len(labels)))
ax2.set_xticklabels(labels, rotation=45, ha='right', fontsize=11, fontstyle='italic')
ax2.set_ylabel('Number of independent\nviral clusters', fontsize=12)
ax2.set_title('Viral dispersal across tree', fontsize=13, fontweight='bold')
ax2.grid(axis='y', alpha=0.3)

# legend
legend_el = [Patch(facecolor=CLASS_COLORS[c], label=c) for c in ['Metabolic','Physiological','Regulatory']]
ax1.legend(handles=legend_el, loc='upper right', fontsize=10, title='AVG class')

plt.tight_layout()
plt.savefig('hgt_cohesion.pdf', bbox_inches='tight', dpi=300)
plt.savefig('hgt_cohesion.png', bbox_inches='tight', dpi=150)
plt.close()
print("saved hgt_cohesion")

# ============ Panel B: donor genus heatmap ============
# collect all genera
all_genera={}
for d in data:
    for g,c in d['donors'].items():
        all_genera[g]=all_genera.get(g,0)+c
top_genera=[g for g,c in sorted(all_genera.items(),key=lambda x:-x[1])[:12]]

# build matrix: gene x genus (fraction of that gene's viruses)
matrix=np.zeros((len(genes_ord), len(top_genera)))
for i,d in enumerate(genes_ord):
    total=sum(d['donors'].values())
    for j,g in enumerate(top_genera):
        matrix[i,j]=100*d['donors'].get(g,0)/total if total>0 else 0

fig, ax = plt.subplots(figsize=(12, 7))
im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd')
ax.set_xticks(range(len(top_genera)))
ax.set_xticklabels(top_genera, rotation=45, ha='right', fontsize=10, fontstyle='italic')
# phylum colour strip (bottom)
phylum_colors={'Firmicutes':'#44BB99','Bacteroidetes':'#EE8866','Spirochaetes':'#BBCC33','Actinobacteria':'#FFAABB','Proteobacteria':'#99DDFF','Other':'#DDDDDD'}
for j,g in enumerate(top_genera):
    ph=phylum_of(g)
    ax.add_patch(plt.Rectangle((j-0.5, -0.5-1.0), 1, 0.5, facecolor=phylum_colors.get(ph,'#DDD'), clip_on=False, edgecolor='white', linewidth=0.5))
ax.set_yticks(range(len(genes_ord)))
ax.set_yticklabels([d['gene'] for d in genes_ord], fontsize=11, fontstyle='italic')
# class colour block to the left of gene name
for i,d in enumerate(genes_ord):
    ax.add_patch(plt.Rectangle((-1.2,i-0.4),0.5,0.8,facecolor=CLASS_COLORS[d['cls']],
                               clip_on=False,edgecolor='none'))
ax.set_xlim(-1.5, len(top_genera)-0.5)
ax.set_ylim(len(genes_ord)-0.5, -2.0)
cbar=plt.colorbar(im, ax=ax, label='% of nearest-neighbor bacteria')
ax.set_title('Donor genera of AVGs', 
             fontsize=13, fontweight='bold', pad=42)
ax.set_xlabel('Bacterial genus (donor candidate)', fontsize=11)

class_leg = [Patch(facecolor=CLASS_COLORS[c], label=c) for c in ['Metabolic','Physiological','Regulatory']]
seen_phyla=[]
for g in top_genera:
    ph=phylum_of(g)
    if ph not in seen_phyla: seen_phyla.append(ph)
phylum_leg=[Patch(facecolor=phylum_colors[p],label=p) for p in seen_phyla]
leg1=ax.legend(handles=class_leg, loc='upper left', bbox_to_anchor=(1.16,1.0), fontsize=9, title='AVG class', frameon=True)
ax.add_artist(leg1)
ax.legend(handles=phylum_leg, loc='upper left', bbox_to_anchor=(1.16,0.62), fontsize=9, title='Phylum', frameon=True)

plt.tight_layout()
plt.savefig('hgt_donors.pdf', bbox_inches='tight', dpi=300)
plt.savefig('hgt_donors.png', bbox_inches='tight', dpi=150)
plt.close()
print("saved hgt_donors")
