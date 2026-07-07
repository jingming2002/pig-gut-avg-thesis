#!/usr/bin/env python3
"""Three-layer Sankey v2: more viral families (top15), coloured by viral type"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.path import Path
import matplotlib.patches as patches

mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

# ===== Data =====
def load_pred(f):
    rows=[]
    with open(f) as fh:
        next(fh)
        for line in fh:
            p=line.rstrip('\n').split(',')
            rows.append({'vOTU':p[0].split('|')[0], 'family':p[2].replace('_like','')})
    return pd.DataFrame(rows)

pg = pd.concat([load_pred(f) for f in [
    'phagcn3_high_amg/final_amg_prediction.csv',
    'phagcn3_high_apg/final_apg_prediction.csv',
    'phagcn3_high_areg/final_areg_prediction.csv']]).drop_duplicates('vOTU')

avg = pd.read_csv('avg_with_function.tsv', sep='\t')
CLS={'metabolic':'AMG','physiological':'APG','regulatory':'AReG'}
avg['avg_class']=avg['Protein Classification'].map(CLS)
vf = avg.groupby('vOTU').agg(
    function=('function', lambda x:x.value_counts().index[0]),
    avg_class=('avg_class', lambda x:x.value_counts().index[0])).reset_index()
m = vf.merge(pg, on='vOTU', how='inner')
mc = m[~m['function'].isin(['(no annotation)','Other'])].copy()
top_funcs = mc['function'].value_counts().head(10).index.tolist()
mc = mc[mc['function'].isin(top_funcs)]

# top 15 families
TOPN_FAM = 15
top_fams = mc['family'].value_counts().head(TOPN_FAM).index.tolist()
mc['family_grp'] = mc['family'].apply(lambda x: x if x in top_fams else 'Other')

# Viral family classification: phage vs eukaryotic virus (suspect)
EUKARYOTIC = {'Geminiviridae','Solemoviridae','Smacoviridae','Caulimoviridae',
              'Virgaviridae','Betaflexiviridae','Secoviridae','Triavirus',
              'Adenoviridae','Papillomaviridae','Paramyxoviridae','Rhabdoviridae',
              'Genomoviridae','Mymonaviridae','Alexandravirus','Fernvirus','Mimiviridae'}

classes = ['AMG','APG','AReG']
functions = mc['function'].value_counts().index.tolist()
families = mc['family_grp'].value_counts().index.tolist()
if 'Other' in families:
    families.remove('Other'); families.append('Other')

flow1 = mc.groupby(['avg_class','function']).size().to_dict()
flow2 = mc.groupby(['function','family_grp']).size().to_dict()

fig, ax = plt.subplots(figsize=(15, 12))
COL_X = {'class':0.05, 'func':0.5, 'fam':0.95}
NODE_W = 0.022; GAP = 0.009
class_colors = {'AMG':'#E07A5F','APG':'#9B2226','AReG':'#3D5A80'}
func_color = '#6A994E'
fam_phage = '#5B8C5A'    # phagefamily green
fam_euk = '#B0A8B9'      # eukaryoticvirus(suspect) grey-purple

def fam_color_of(fam):
    if fam == 'Other': return '#CCCCCC'
    return fam_euk if fam in EUKARYOTIC else fam_phage

def layout_column(items, sizes):
    total = sum(sizes[it] for it in items)
    h_total = 1.0 - GAP*(len(items)+1)
    pos = {}; y = GAP
    for it in items:
        h = h_total * sizes[it]/total
        pos[it] = (y, y+h); y += h + GAP
    return pos

class_sizes = mc['avg_class'].value_counts().to_dict()
func_sizes = mc['function'].value_counts().to_dict()
fam_sizes = mc['family_grp'].value_counts().to_dict()
class_pos = layout_column(classes, class_sizes)
func_pos = layout_column(functions, func_sizes)
fam_pos = layout_column(families, fam_sizes)

def draw_node(x, ypos, color, label, align='left', fs=8):
    yb, yt = ypos
    ax.add_patch(patches.Rectangle((x,yb),NODE_W,yt-yb,facecolor=color,edgecolor='white',linewidth=0.5,zorder=3))
    ymid=(yb+yt)/2
    if align=='left':
        ax.text(x-0.007,ymid,label,ha='right',va='center',fontsize=fs,zorder=4)
    else:
        ax.text(x+NODE_W+0.007,ymid,label,ha='left',va='center',fontsize=fs,zorder=4)

for c in classes:
    draw_node(COL_X['class'], class_pos[c], class_colors[c], f"{c} ({class_sizes[c]})", 'left', 10)
for f in functions:
    yb,yt=func_pos[f]; ymid=(yb+yt)/2
    ax.add_patch(patches.Rectangle((COL_X['func'],yb),NODE_W,yt-yb,facecolor=func_color,edgecolor='white',linewidth=0.5,zorder=3))
    ax.text(COL_X['func']+NODE_W+0.007,ymid,f"{f} ({func_sizes[f]})",ha='left',va='center',fontsize=8.5,zorder=4)
for fam in families:
    draw_node(COL_X['fam'], fam_pos[fam], fam_color_of(fam), f"{fam} ({fam_sizes[fam]})", 'right', 8)

def draw_flow(x0,y0r,x1,y1r,color,alpha=0.3):
    y0b,y0t=y0r; y1b,y1t=y1r; xc=(x0+x1)/2
    verts=[(x0,y0b),(xc,y0b),(xc,y1b),(x1,y1b),(x1,y1t),(xc,y1t),(xc,y0t),(x0,y0t),(x0,y0b)]
    codes=[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY]
    ax.add_patch(patches.PathPatch(Path(verts,codes),facecolor=color,edgecolor='none',alpha=alpha,zorder=1))

class_used={c:class_pos[c][0] for c in classes}
func_used_left={f:func_pos[f][0] for f in functions}
func_used_right={f:func_pos[f][0] for f in functions}
fam_used={fam:fam_pos[fam][0] for fam in families}

for c in classes:
    for f in functions:
        cnt=flow1.get((c,f),0)
        if cnt<3: continue
        cb,ct=class_pos[c]; h_c=(ct-cb)*cnt/class_sizes[c]
        fb,ft=func_pos[f]; h_f=(ft-fb)*cnt/func_sizes[f]
        draw_flow(COL_X['class']+NODE_W,(class_used[c],class_used[c]+h_c),COL_X['func'],(func_used_left[f],func_used_left[f]+h_f),class_colors[c],0.3)
        class_used[c]+=h_c; func_used_left[f]+=h_f

for f in functions:
    for fam in families:
        cnt=flow2.get((f,fam),0)
        if cnt<3: continue
        fb,ft=func_pos[f]; h_f=(ft-fb)*cnt/func_sizes[f]
        mb,mt=fam_pos[fam]; h_m=(mt-mb)*cnt/fam_sizes[fam]
        c = fam_euk if fam in EUKARYOTIC else func_color
        draw_flow(COL_X['func']+NODE_W,(func_used_right[f],func_used_right[f]+h_f),COL_X['fam'],(fam_used[fam],fam_used[fam]+h_m),c,0.2)
        func_used_right[f]+=h_f; fam_used[fam]+=h_m

ax.text(COL_X['class']+NODE_W/2,1.03,'AVG class',ha='center',fontsize=12,fontweight='bold')
ax.text(COL_X['func']+NODE_W/2,1.03,'Function',ha='center',fontsize=12,fontweight='bold')
ax.text(COL_X['fam']+NODE_W/2,1.03,'Viral family',ha='center',fontsize=12,fontweight='bold')

# legend: phage vs eukaryotic virus
leg = [patches.Patch(facecolor=fam_phage, label='Bacteriophage family'),
       patches.Patch(facecolor=fam_euk, label='Eukaryotic virus (likely misclassified)')]
ax.legend(handles=leg, loc='lower center', bbox_to_anchor=(0.5,-0.06), ncol=2, fontsize=9, frameon=False)

ax.set_xlim(-0.2,1.2); ax.set_ylim(-0.04,1.07); ax.axis('off')
ax.set_title(f"AVG class -> Function -> Viral family (n={len(mc)} vOTUs, pig gut viruses)",fontsize=14,fontweight='bold',pad=28)
plt.savefig('avg_sankey_v2.pdf', bbox_inches='tight', dpi=300)
plt.savefig('avg_sankey_v2.png', bbox_inches='tight', dpi=200)
print(f"saved, n={len(mc)}, families shown={len(families)}")
