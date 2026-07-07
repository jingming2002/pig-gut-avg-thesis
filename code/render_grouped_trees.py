#!/usr/bin/env python3
"""Render trees grouped by AVG class: 3 figures (metabolic/physiological/regulatory), each with all trees side by side"""
from Bio import Phylo
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
}
PHYLUM_COLORS={'Firmicutes':'#2CA25F','Bacteroidetes':'#E6550D','Spirochaetes':'#BCBD22',
               'Actinobacteria':'#E377C2','Other':'#888888'}
VIRUS_COLOR='#D62728'

def leaf_color(name):
    if name.startswith('VIRUS_'): return VIRUS_COLOR
    if name.startswith('BACT_'):
        return PHYLUM_COLORS.get(PHYLUM.get(name.split('_')[1],'Other'),'#888888')
    return '#000000'

def draw_one(ax, treepath, subtitle):
    tree = Phylo.read(treepath, "newick")
    n = tree.count_terminals()
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: None, show_confidence=False)
    leaves = tree.get_terminals()
    depths = tree.depths()
    for i, leaf in enumerate(leaves):
        x = depths[leaf]; y = i + 1
        ax.scatter([x],[y], s=10, c=leaf_color(leaf.name), zorder=5, edgecolors='none')
    ax.set_title(subtitle, fontsize=12, fontweight='bold')
    ax.axis('off')
    return n

# three-class config
GROUPS = {
    'Metabolic': [
        ("cysH_sulfur_metabolism_K00390","cysH (Sulfur)"),
        ("phoH_phosphate_metabolism_K06217","phoH (Phosphate)"),
        ("folD_folate_C1_metabolism_K01491","folD (One-carbon)"),
        ("iscS_cysteine_desulfurase_K04487","iscS (Fe-S)"),
    ],
    'Physiological': [
        ("mazF_toxin_antitoxin_K07171","mazF (Toxin-AT)"),
        ("spoVT_sporulation_regulation_K04769","spoVT (Sporulation)"),
    ],
    'Regulatory': [
        ("lexA_SOS_lysogeny_repressor_K01356","lexA (SOS/lysogeny)"),
        ("dicA_regulatory_K22300","dicA (Cell-div reg)"),
        ("sinR_regulatory_K19449","sinR (Biofilm reg)"),
    ],
}

legend_el=[
    Line2D([0],[0],marker='o',color='w',markerfacecolor=VIRUS_COLOR,markersize=10,label='Virus (AVG)'),
    Line2D([0],[0],marker='o',color='w',markerfacecolor=PHYLUM_COLORS['Firmicutes'],markersize=10,label='Firmicutes'),
    Line2D([0],[0],marker='o',color='w',markerfacecolor=PHYLUM_COLORS['Bacteroidetes'],markersize=10,label='Bacteroidetes'),
    Line2D([0],[0],marker='o',color='w',markerfacecolor=PHYLUM_COLORS['Spirochaetes'],markersize=10,label='Spirochaetes'),
    Line2D([0],[0],marker='o',color='w',markerfacecolor=PHYLUM_COLORS['Actinobacteria'],markersize=10,label='Actinobacteria'),
]

for cls, gene_list in GROUPS.items():
    ntrees = len(gene_list)
    fig, axes = plt.subplots(1, ntrees, figsize=(5.5*ntrees, 11))
    if ntrees == 1: axes = [axes]
    
    for ax, (gene, subtitle) in zip(axes, gene_list):
        path = f"phagcn3_high_apg/trees_v2/{gene}/tree.treefile"
        try:
            n = draw_one(ax, path, subtitle)
        except Exception as e:
            ax.text(0.5,0.5,f"failed:\n{e}",ha='center',transform=ax.transAxes)
            ax.axis('off')
    
    fig.suptitle(f"{cls} AVGs — viral phylogenetic placement", fontsize=15, fontweight='bold', y=1.00)
    # shared legend at bottom
    fig.legend(handles=legend_el, loc='lower center', ncol=5, fontsize=10, frameon=True, bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout()
    out = f"trees_{cls.lower()}.png"
    plt.savefig(out, dpi=130, bbox_inches='tight')
    plt.savefig(f"trees_{cls.lower()}.pdf", bbox_inches='tight')
    plt.close()
    print(f"saved {out} ({ntrees} trees)")

print("done - 3 grouped figures")
