# Auxiliary Viral Genes in the Pig Gut Virome

MSc thesis project — Wang Jingming
Supervisors: Dr Anne Kupczok, Dr Liu Jun (Wageningen University)

## Overview

This project investigates auxiliary viral genes (AVGs) in the pig gut virome —
virus-encoded genes that augment host metabolism (AMG), physiology (APG), or
regulation (AReG). The central question is whether the three AVG classes represent
three dimensions of viral adaptation, examined through their function, horizontal
gene transfer (HGT) history, and predicted bacterial hosts.

## Key findings

- Three AVG classes map to three functional dimensions (nutrient acquisition,
  defense, host manipulation); function is predicted by class, not viral taxonomy.
- AVGs are acquired predominantly by HGT from Firmicutes donors, corroborated by
  host prediction (hosts are 65-74% Firmicutes) - donor approximates host.
- A cohesion gradient: metabolic genes are ancestral/clustered, regulatory genes
  are dispersed via recurrent recent HGT.

## Pipeline

PVD viral contigs -> CheckV (completeness >90%) -> CheckAMG (AVG identification)
-> three branches: PhaGCN3 (taxonomy), iPHoP (host prediction), DIAMOND blastp
(vs clean bacterial refs) -> MAFFT -> trimAl -> IQ-TREE -> ggtree (HGT detection)

## Repository structure

- code/               Analysis scripts
- results/trees/       Final ML trees for the 9 AVG genes (.treefile)
- results/avg_predictions/   PhaGCN3 high-confidence AVG predictions
- DATA_SOURCES.txt     Where to obtain the public databases used
- *_environment.yml    Conda environments

## The nine AVG genes

Metabolic:     cysH (K00390), phoH (K06217), folD (K01491), iscS (K04487)
Physiological: mazF (K07171), spoVT (K04769)
Regulatory:    lexA (K01356), dicA (K22300), sinR (K19449)

## Key parameters

- CheckV: completeness > 90%
- CheckAMG: min-flank-Vscore >= 10, contig >= 5 kb, >= 4 ORFs
- PhaGCN3: high-confidence predictions only
- DIAMOND: e-value <= 1e-5, max-target-seqs 50
- CD-HIT: c = 0.9 (dicA: c = 0.8)
- IQ-TREE: -m LG+G -bb 1000 -alrt 1000
- HGT support: SH-aLRT >= 80 & UFBoot >= 95

## Data

Large public databases are not included; see DATA_SOURCES.txt for how to obtain
them. Large intermediate files remain on the Wageningen HPC and can be regenerated
with the scripts above.
