#!/usr/bin/env Rscript
# Circular tree + bootstrap node symbols (black dots UFBoot>=95) + dicA downsampling
suppressMessages({
  library(ggtree); library(ggplot2); library(treeio); library(patchwork); library(ape)
})

phylum_map <- c(
  Clostridium="Firmicutes",Ruminococcus="Firmicutes",Blautia="Firmicutes",
  Faecalibacterium="Firmicutes",Streptococcus="Firmicutes",Lactobacillus="Firmicutes",
  Roseburia="Firmicutes",Oscillibacter="Firmicutes",Agathobacter="Firmicutes",
  Megasphaera="Firmicutes",Phascolarctobacterium="Firmicutes",Eubacterium="Firmicutes",
  Coprococcus="Firmicutes",Dorea="Firmicutes",Lachnoclostridium="Firmicutes",
  Anaerostipes="Firmicutes",Subdoligranulum="Firmicutes",
  Bacteroides="Bacteroidetes",Prevotella="Bacteroidetes",Parabacteroides="Bacteroidetes",
  Alistipes="Bacteroidetes",Treponema="Spirochaetes",
  Bifidobacterium="Actinobacteria",Collinsella="Actinobacteria"
)
phylum_colors <- c(
  "Virus"="#D62728","Firmicutes"="#2CA25F","Bacteroidetes"="#E6550D",
  "Spirochaetes"="#BCBD22","Actinobacteria"="#E377C2","Other"="#888888"
)

classify_tip <- function(name) {
  if (startsWith(name,"VIRUS_")) return("Virus")
  if (startsWith(name,"BACT_")) {
    genus <- strsplit(name,"_")[[1]][2]
    ph <- phylum_map[genus]
    if (is.na(ph)) return("Other")
    return(unname(ph))
  }
  return("Other")
}

downsample_tree <- function(tr, max_total=250) {
  if (length(tr$tip.label) <= max_total) return(tr)
  tips <- tr$tip.label
  v_tips <- tips[startsWith(tips,"VIRUS_")]
  b_tips <- tips[startsWith(tips,"BACT_")]
  frac <- max_total/length(tips)
  set.seed(42)
  keep <- c(sample(v_tips, max(1,round(length(v_tips)*frac))),
            sample(b_tips, max(1,round(length(b_tips)*frac))))
  drop.tip(tr, setdiff(tips, keep))
}

make_circ <- function(treepath, subtitle, downsample=FALSE) {
  tr <- read.tree(treepath)
  norig <- length(tr$tip.label)
  if (downsample) {
    tr <- downsample_tree(tr, 250)
    subtitle <- paste0(subtitle, "\n(", length(tr$tip.label), " of ", norig, " tips)")
  }
  
  # parse node label UFBoot (format aLRT/UFBoot)
  # tr$node.label e.g. "94.8/100"
  ufboot <- rep(NA, length(tr$node.label))
  for (i in seq_along(tr$node.label)) {
    lab <- tr$node.label[i]
    if (!is.na(lab) && grepl("/", lab)) {
      parts <- strsplit(lab, "/")[[1]]
      ufboot[i] <- suppressWarnings(as.numeric(parts[2]))
    }
  }
  
  tip_class <- sapply(tr$tip.label, classify_tip)
  dat <- data.frame(label=tr$tip.label, group=tip_class)
  
  p <- ggtree(tr, layout="fan", open.angle=12, size=0.2) %<+% dat +
    geom_tippoint(aes(color=group), size=1.3) +
    scale_color_manual(values=phylum_colors, name="") +
    ggtitle(subtitle) +
    theme(plot.title=element_text(size=12,face="bold",hjust=0.5),
          legend.position="none")
  
  # add bootstrap node symbols: black dots for UFBoot>=95
  # ggtree internal node numbering = Ntip+1 .. Ntip+Nnode
  ntip <- length(tr$tip.label)
  nodedf <- data.frame(node=(ntip+1):(ntip+tr$Nnode), ufboot=ufboot)
  nodedf_strong <- nodedf[!is.na(nodedf$ufboot) & nodedf$ufboot>=95, ]
  
  if (nrow(nodedf_strong) > 0) {
    p <- p + geom_point2(aes(subset=(node %in% nodedf_strong$node)),
                         size=0.8, color="black", shape=16)
  }
  return(p)
}

groups <- list(
  Metabolic=list(
    list("cysH_sulfur_metabolism_K00390","cysH (Sulfur)",FALSE),
    list("phoH_phosphate_metabolism_K06217","phoH (Phosphate)",FALSE),
    list("folD_folate_C1_metabolism_K01491","folD (One-carbon)",FALSE),
    list("iscS_cysteine_desulfurase_K04487","iscS (Fe-S)",FALSE)
  ),
  Physiological=list(
    list("mazF_toxin_antitoxin_K07171","mazF (Toxin-AT)",FALSE),
    list("spoVT_sporulation_regulation_K04769","spoVT (Sporulation)",FALSE)
  ),
  Regulatory=list(
    list("lexA_SOS_lysogeny_repressor_K01356","lexA (SOS/lysogeny)",FALSE),
    list("dicA_regulatory_K22300","dicA (Cell-div reg)",TRUE),
    list("sinR_regulatory_K19449","sinR (Biofilm reg)",FALSE)
  )
)

for (cls in names(groups)) {
  gl <- groups[[cls]]
  plots <- lapply(gl, function(item){
    make_circ(paste0("phagcn3_high_apg/trees_v2/",item[[1]],"/tree.treefile"),
              item[[2]], item[[3]])
  })
  ncol_use <- length(plots)
  combined <- wrap_plots(plots, ncol=ncol_use) +
    plot_annotation(title=paste0(cls," AVGs — viral phylogenetic placement (black dots: UFBoot>=95)"),
                    theme=theme(plot.title=element_text(size=14,face="bold",hjust=0.5)))
  out <- paste0("circboot_",tolower(cls),".png")
  ggsave(out, combined, width=5.5*ncol_use, height=6.5, dpi=140)
  cat("saved",out,"\n")
}
cat("done\n")
