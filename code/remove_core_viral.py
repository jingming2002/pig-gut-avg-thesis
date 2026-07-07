#!/usr/bin/env python3
"""
removethree classes AVG clearly core invirusgene(based on Pfam Name),
countcleandistribution

basis: Martin et al. 2025 - AVG defined as "auxiliary to core viral functions",
core replication/structural/packaginggene
"""

import pandas as pd
import re

df = pd.read_csv("checkamg_output/results/final_results.tsv", sep="\t", low_memory=False)
df_high = df[df['Protein Viral Origin Confidence'] == 'high'].copy()

# ============================================================
# core viral gene Pfam blacklist (exact Pfam Name match)
# remove only genes clearly viral replication/structural/packaging/assembly
# ============================================================
CORE_VIRAL_PFAMS = [
    # ssDNA-binding / replication (SSB type)
    "Single-strand binding protein family",
    "Phage ORF5 protein",
    "Putative phage ssDNA-binding domain",
    "Enterobacter phage Enc34, ssDNA-binding protein",
    "Putative DnaT-like ssDNA binding protein",
    "Single-stranded DNA-binding protein, Bacteriophage T7",
    "gp32 DNA binding protein like",
    "Phage Single-stranded DNA-binding protein",
    "Viral DNA-binding protein, zinc binding domain",
    # replicative helicase/recombination
    "SNF2-related domain",
    # structural assembly chaperone
    "Chaperone of endosialidase",
]

# also use keyword fallback to catch core structural/packaging genes missed in other classes
CORE_VIRAL_KEYWORDS = [
    r'terminase', r'capsid', r'\bportal\b', r'baseplate', r'tail fiber',
    r'tail spike', r'tape.?measure', r'major.*coat', r'prohead', r'tail tube',
    r'tail terminator', r'tail assembly', r'head.?tail', r'\bsheath\b',
    r'holin', r'endolysin', r'\bspanin\b', r'portal protein',
]
core_kw_pattern = '|'.join(CORE_VIRAL_KEYWORDS)

def is_core_viral(pfam_name):
    if pd.isna(pfam_name):
        return False
    name = str(pfam_name)
    # exact blacklist
    if name in CORE_VIRAL_PFAMS:
        return True
    # keyword fallback
    if re.search(core_kw_pattern, name.lower()):
        return True
    return False

# flag core viral genes
df_high['is_core_viral'] = df_high['Pfam Name'].apply(is_core_viral)

# ============================================================
# summarise removals
# ============================================================
print("="*80)
print("corevirusgeneremovecount")
print("="*80)
for cls in ['metabolic', 'physiological', 'regulatory']:
    sub = df_high[df_high['Protein Classification'] == cls]
    core = sub['is_core_viral'].sum()
    print(f"  {cls:15s}: remove {core:5d} / {len(sub):5d} ({100*core/len(sub):4.1f}%)  → keep {len(sub)-core}")

total_core = df_high['is_core_viral'].sum()
print(f"\n  totalremove: {total_core} / {len(df_high)} ({100*total_core/len(df_high):.1f}%)")

# specific removed Pfams (verification)
print("\nwasremovecore ofvirusgene (by Pfam):")
removed = df_high[df_high['is_core_viral']]
for name, c in removed['Pfam Name'].value_counts().items():
    print(f"  {c:5d}  {str(name)[:60]}")

# ============================================================
# after removal, recompute functional classification
# ============================================================
df_clean = df_high[~df_high['is_core_viral']].copy()

# full functional classification rules (as before)
CATEGORY_RULES = [
    ("Anti-defense / Anti-restriction", [r'anti.?restriction', r'ardA', r'ardC', r'anti.?defense', r'thoeris anti', r'anti.?crispr', r'acrI', r'dna mimic', r'\bocr\b', r'darb']),
    ("Defense system (RM/CRISPR/Abi)", [r'restriction', r'crispr', r'\bcas\b', r'\babi', r'abortive', r'thoeris', r'defense', r'retron', r'muts']),
    ("Toxin-antitoxin", [r'toxin', r'antitoxin', r'mazf', r'maze', r'pemk', r'hica', r'hicb', r'hipa', r'hipb', r'rele', r'relb', r'pare', r'pard', r'vapc', r'vapb', r'abieii', r'\bdoc\b', r'\bphd\b', r'pemi']),
    ("Transcriptional regulation (HTH)", [r'helix.?turn.?helix', r'\bhth\b', r'cro/c1', r'cro\b', r'winged helix', r'transcriptional regulat', r'sigma.?70', r'sigma factor', r'repressor', r'antirepressor', r'anti.?repressor', r'\bbro\b', r'kilac', r'arac', r'tetr', r'lysr', r'merr', r'xre', r'lexa', r'\bant\b']),
    ("Peptidase / Protease (regulatory)", [r'peptidase s24', r'sos response', r'peptidase', r'protease']),
    ("DNA-binding (host-like)", [r'bacterial dna.?binding', r'\bhu\b', r'\bihf\b', r'parb', r'\bpara\b', r'histone']),
    ("DNA replication / repair", [r'helicase', r'primase', r'polymerase', r'recombinase', r'resolvase', r'topoisomerase', r'nuclease', r'replication', r'holliday', r'\bdut']),
    ("DNA methylation / modification", [r'methyltransferase', r'methylase', r'\bdam\b', r'\bdcm\b', r'n-6 dna', r'c-5 cytosine']),
    ("Integration / Transposition", [r'integrase', r'transposase', r'invertase', r'\bxer', r'recombinase']),
    ("Chaperone / Heat shock", [r'chaperon', r'groel', r'groes', r'cpn60', r'cpn10', r'\bhsp', r'tcp-1']),
    ("Sulfur metabolism", [r'phosphoadenosine phosphosulfate', r'\bpaps\b', r'sulf', r'cysh', r'cysteine synth', r'rubrerythrin']),
    ("Phosphate metabolism", [r'phoh', r'phosphate', r'pyrophosphatase', r'phosphatase', r'pts hpr']),
    ("Nucleotide / Folate metabolism", [r'ribonucleotide reductase', r'\bnrd', r'thymidylate', r'tetrahydrofolate', r'folate', r'nucleotide', r'\bdut', r'purine', r'pyrimidine']),
    ("Amino acid metabolism", [r'aminotransferase', r'glutam', r'gamma-glutamyl', r'aig2', r'amino acid']),
    ("Acyltransferase / Acetyltransferase", [r'acetyltransferase', r'gnat', r'hexapeptide', r'acyltransferase', r'gdsl']),
    ("Cofactor / Vitamin", [r'cobalamin', r'cobt', r'biotin', r'thiamine', r'riboflavin']),
    ("Oxidoreductase / Electron transfer", [r'oxidoreductase', r'reductase', r'dehydrogenase', r'oxygenase', r'oxidase', r'ferredoxin', r'ferritin', r'radical sam', r'\bnad']),
    ("Transport", [r'abc transporter', r'transporter', r'permease', r'\bompa\b', r'porin']),
    ("Sporulation", [r'sporulation', r'\bspov', r'\bspo0', r'spore']),
    ("Motility", [r'flagell', r'\bfla', r'\bflg', r'\bfli', r'chemotaxis', r'pilus']),
    ("Cell wall / Membrane", [r'peptidoglycan', r'murein', r'\bmur', r'cell wall', r'lysozyme', r'amidase', r'sortase']),
    ("Signal transduction", [r'kinase', r'two.?component', r'response regulator', r'histidine kinase', r'\bfic\b']),
    ("Secretion system", [r'secretion', r'vgrg', r't6ss', r'type vi', r'spike']),
    ("Virulence-associated", [r'virulence', r'vape', r'large polyvalent']),
]

def classify(text):
    if pd.isna(text) or str(text).strip()=="":
        return "(no Pfam)"
    t = str(text).lower()
    for cat, pats in CATEGORY_RULES:
        for p in pats:
            if re.search(p, t):
                return cat
    return "Other"

print("\n" + "="*80)
print("removecorevirusafter genes - three classescleandistribution")
print("="*80)

results = {}
for cls in ['metabolic', 'physiological', 'regulatory']:
    sub = df_clean[df_clean['Protein Classification'] == cls].copy()
    sub['cat'] = sub['Pfam Name'].apply(classify)
    cc = sub['cat'].value_counts()
    results[cls] = cc
    print(f"\n【{cls.upper()}】removeafter {len(sub)} ")
    for cat, c in cc.head(12).items():
        print(f"  {c:5d} ({100*c/len(sub):4.1f}%)  {cat}")

# save
summary = pd.DataFrame(results).fillna(0).astype(int)
summary.to_csv("avg_function_distribution_clean.tsv", sep="\t")
df_clean.to_csv("avg_high_conf_clean.tsv", sep="\t", index=False)
print(f"\nalreadysave: avg_function_distribution_clean.tsv (distribution)")
print(f"alreadysave: avg_high_conf_clean.tsv (removecoreafter genesclean AVG set)")
EOF
