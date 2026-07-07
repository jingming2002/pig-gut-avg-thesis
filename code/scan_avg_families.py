#!/usr/bin/env python3
"""
scan high-confidence AVG annotation of,countgenedistribution
treegene

:
  - CheckAMG final_results.tsv or similar annotation table
  - high-confidence AVG protein FASTA (optional, for length stats)

:
  - KEGG KO abundance table (top 50)
  - Pfam domain abundance table (top 50)
  - functional category enrichment (KEGG pathway)
  - Candidate gene family recommendations(by x )
"""

import pandas as pd
import argparse
from collections import Counter
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser(description="Scan AVG annotation and recommend gene families for phylogenetic analysis")
    p.add_argument("--annotation", required=True, help="CheckAMG final_results.tsv or gene_annotations.tsv")
    p.add_argument("--output-prefix", default="avg_family_scan", help="Output file prefix")
    p.add_argument("--top-n", type=int, default=50, help="Show top N families by abundance")
    return p.parse_args()

def main():
    args = parse_args()
    
    # read annotation table
    print("[1/5] read...")
    try:
        # try reading final_results.tsv (CheckAMG output)
        df = pd.read_csv(args.annotation, sep="\t", low_memory=False)
    except:
        # if failed, try gene_annotations.tsv
        df = pd.read_csv(args.annotation, sep=",", low_memory=False)
    
    print(f"      total: {len(df):,}")
    print(f"      column: {df.columns.tolist()}")
    
    # ---- 1. KEGG KO counts ----
    print("\n[2/5] count KEGG KO ...")
    kegg_cols = [c for c in df.columns if 'kegg' in c.lower() or 'ko' in c.lower()]
    
    if kegg_cols:
        kegg_col = kegg_cols[0]
        print(f"      usecolumn: {kegg_col}")
        kegg_counts = df[kegg_col].fillna("").str.split(";").explode().value_counts()
        kegg_counts = kegg_counts[kegg_counts.index != ""]
        
        kegg_df = pd.DataFrame({
            "KEGG_KO": kegg_counts.index[:args.top_n],
            "count": kegg_counts.values[:args.top_n]
        })
        kegg_df.to_csv(f"{args.output_prefix}_kegg_top{args.top_n}.tsv", sep="\t", index=False)
        print(f"      Top {args.top_n} KO: alreadysave")
        print(kegg_df.head(10))
    else:
        print("      not found KEGG column,skip")
        kegg_df = None
    
    # ---- 2. Pfam domain counts ----
    print("\n[3/5] count Pfam domain ...")
    pfam_cols = [c for c in df.columns if 'pfam' in c.lower()]
    
    if pfam_cols:
        pfam_col = pfam_cols[0]
        print(f"      usecolumn: {pfam_col}")
        pfam_counts = df[pfam_col].fillna("").str.split(";").explode().value_counts()
        pfam_counts = pfam_counts[pfam_counts.index != ""]
        
        pfam_df = pd.DataFrame({
            "Pfam_domain": pfam_counts.index[:args.top_n],
            "count": pfam_counts.values[:args.top_n]
        })
        pfam_df.to_csv(f"{args.output_prefix}_pfam_top{args.top_n}.tsv", sep="\t", index=False)
        print(f"      Top {args.top_n} Pfam: alreadysave")
        print(pfam_df.head(10))
    else:
        print("      not found Pfam column,skip")
        pfam_df = None
    
    # ---- 3. AVG class distribution ----
    print("\n[4/5] AVG class distribution...")
    class_cols = [c for c in df.columns if 'classification' in c.lower() or 'type' in c.lower()]
    if class_cols:
        class_col = class_cols[0]
        class_counts = df[class_col].value_counts()
        print(class_counts)
        class_counts.to_csv(f"{args.output_prefix}_classification.tsv", sep="\t", header=["count"])
    
    # ---- 4. Candidate gene family recommendations ----
    print("\n[5/5] gene...")
    print("\n" + "="*80)
    print(" 5-10 tree buildinggene:")
    print("="*80)
    
    candidates = []
    
    # high-abundance families from KEGG
    if kegg_df is not None:
        candidates_kegg = kegg_df[kegg_df["count"] >= 20].head(5)
        print("\nfrom KEGG KO (sample size >= 20)")
        for idx, row in candidates_kegg.iterrows():
            print(f"  - {row['KEGG_KO']}: {row['count']} ")
            candidates.append(f"KEGG:{row['KEGG_KO']}")
    
    # high-abundance families from Pfam
    if pfam_df is not None:
        candidates_pfam = pfam_df[pfam_df["count"] >= 30].head(5)
        print("\nfrom Pfam domain (sample size >= 30)")
        for idx, row in candidates_pfam.iterrows():
            print(f"  - {row['Pfam_domain']}: {row['count']} ")
            candidates.append(f"Pfam:{row['Pfam_domain']}")
    
    print("\n(already)")
    manual_candidates = [
        ("K02588", "nifH (nitrogen fixation)"),
        ("K04561", "nirK (denitrification)"),
        ("K00266", "glnA (glutamine synthetase)"),
        ("PF00004", "AAA ATPase domain"),
        ("PF00023", "Ankyrin repeat"),
    ]
    for ko, desc in manual_candidates:
        if kegg_df is not None and ko in kegg_df["KEGG_KO"].values:
            count = kegg_df[kegg_df["KEGG_KO"] == ko]["count"].values[0]
            print(f"  - {ko}: {desc} ({count} )")
            candidates.append(ko)
    
    # save candidate list
    with open(f"{args.output_prefix}_candidates.txt", "w") as f:
        f.write("\n".join(set(candidates)))
    
    print("\n" + "="*80)
    print("next step:")
    print("  1. fromon 5-10 gene")
    print("  2. eachextract KEGG KO / Pfam ID")
    print("  3. fromprotein FASTA extractsequence")
    print("  4. BLAST/DIAMOND  NR, top hits")
    print("  5. use MAFFT/trimAl/IQ-TREE2 tree building")
    print("="*80)

if __name__ == "__main__":
    main()
