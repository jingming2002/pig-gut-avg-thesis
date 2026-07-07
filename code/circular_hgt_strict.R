#!/usr/bin/env Rscript
# Circular tree: mark only terminal branches of truly nested HGT viruses (no parent-node chaining)
suppressMessages({
  library(ggtree); library(ggplot2); library(treeio); library(patchwork); library(ape)
})

phylum_map <- c(
  Clostridium="Firmicutes",Ruminococcus="Firmicutes",Blautia="Firmicutes",
  Faecalibacterium="Firmicutes",Streptococcus="Firmicutes",Lactobacillus="Firmicutes",
  Roseburia="Firmicutes",Oscillibacter="Firmicutes",Agathobacter="Firmicutes",
  Megasphaera="Firmicutes",Phascolarctobacterium="Firmicutes",Eubacterium="Firmicutes",
  Bacteroides="Bacteroidetes",Prevotella="Bacteroidetes",Parabacteroides="Bacteroidetes",
  Alistipes="Bacteroidetes",Phocaeicola="Bacteroidetes",Segatella="Bacteroidetes",
  Treponema="Spirochaetes",Bifidobacterium="Actinobacteria",Escherichia="Proteobacteria"
)
tip_colors <- c(
  "Firmicutes"="#2CA25F","Bacteroidetes"="#E6550D","Spirochaetes"="#BCBD22",
  "Actinobacteria"="#E377C2","Proteobacteria"="#3182BD","Virus"="#4A4A4A","Other"="#888888"
)

classify_tip <- function(name) {
  if (startsWith(name,"VIRUS_")) return("Virus")
  if (startsWith(name,"BACT_")) {
    g <- strsplit(name,"_")[[1]][2]; ph <- phylum_map[g]
    if (is.na(ph)) return("Other"); return(unname(ph))
  }
  return("Other")
}

downsample_tree <- function(tr, hgt_set, max_total=280) {
  if (length(tr$tip.label) <= max_total) return(tr)
  tips <- tr$tip.label
  v <- tips[startsWith(tips,"VIRUS_")]; b <- tips[startsWith(tips,"BACT_")]
  hgt_v <- intersect(v, hgt_set); other_v <- setdiff(v, hgt_set)
  frac <- max_total/length(tips); set.seed(42)
  keep <- c(hgt_v, sample(other_v,max(1,round(length(other_v)*frac))),
            sample(b,max(1,round(length(b)*frac))))
  drop.tip(tr, setdiff(tips, keep))
}

make_circ <- function(treepath, hgtfile, subtitle, downsample=FALSE) {
  tr <- read.tree(treepath)
  hgt_set <- if (file.exists(hgtfile)) { x<-readLines(hgtfile); x[nchar(x)>0] } else character(0)
  norig <- length(tr$tip.label)
  if (downsample) {
    tr <- downsample_tree(tr, hgt_set, 280)
    subtitle <- paste0(subtitle,"\n(",length(tr$tip.label),"/",norig," tips)")
  }
  ntip <- length(tr$tip.label)
  hgt_in <- intersect(tr$tip.label, hgt_set)

  # ---- Mark only terminal branches of truly nested HGT viruses; never chain-colour parent nodes ----
  n_total <- ntip + tr$Nnode
  edge_group <- rep("normal", n_total)
  for (tipname in hgt_in) {
    idx <- which(tr$tip.label==tipname)
    if (length(idx)) edge_group[idx] <- "HGT"   # mark onlyleavesits own terminal branch
    # do NOT mark parent node! (chain-colouring removed)
  }

  # support values(dual-support black dots)
  alrt <- rep(NA,length(tr$node.label)); ufboot <- rep(NA,length(tr$node.label))
  for (i in seq_along(tr$node.label)) {
    lab <- tr$node.label[i]
    if (!is.na(lab) && grepl("/",lab)) {
      pp <- strsplit(lab,"/")[[1]]
      alrt[i] <- suppressWarnings(as.numeric(pp[1]))
      ufboot[i] <- suppressWarnings(as.numeric(pp[2]))
    }
  }

  p <- ggtree(tr, layout="fan", open.angle=12, size=0.4)
  p$data$edgegrp <- edge_group[p$data$node]
  p <- p + aes(color=edgegrp) +
    scale_color_manual(values=c("HGT"="#E8000B","normal"="#AAAAAA"), guide="none")
  # tip points (by phylum)
  dat <- data.frame(node=1:ntip, label=tr$tip.label)
  dat$tipclass <- sapply(dat$label, classify_tip)
  p <- p + geom_tippoint(data=function(d){
             m <- merge(d, dat, by="node", all.x=TRUE)
             m[!is.na(m$tipclass),]
           }, aes(fill=tipclass), color="black", shape=21, size=1.4, stroke=0.15) +
    scale_fill_manual(values=tip_colors, name="", na.value="transparent")
  # dual-support black dots
  strong <- (ntip+1):(ntip+tr$Nnode)
  strong <- strong[!is.na(ufboot) & !is.na(alrt) & alrt>=80 & ufboot>=95]
  if (length(strong)>0) {
    p <- p + geom_point2(aes(subset=(node %in% strong)), size=0.5, color="black", shape=16)
  }
  n_hgt <- length(hgt_in)
  p <- p + ggtitle(paste0(subtitle, "  [", n_hgt, " HGT]")) +
    theme(plot.title=element_text(size=11.5,face="bold",hjust=0.5), legend.position="none")
  return(p)
}

groups <- list(
  Metabolic=list(
    list("cysH_sulfur_metabolism_K00390","cysH (Sulfur)",FALSE),
    list("phoH_phosphate_metabolism_K06217","phoH (Phosphate)",FALSE),
    list("folD_folate_C1_metabolism_K01491","folD (One-carbon)",FALSE),
    list("iscS_cysteine_desulfurase_K04487","iscS (Fe-S)",FALSE)),
  Physiological=list(
    list("mazF_toxin_antitoxin_K07171","mazF (Toxin-AT)",FALSE),
    list("spoVT_sporulation_regulation_K04769","spoVT (Sporulation)",FALSE)),
  Regulatory=list(
    list("lexA_SOS_lysogeny_repressor_K01356","lexA (SOS/lysogeny)",FALSE),
    list("dicA_regulatory_K22300","dicA (Cell-div reg)",TRUE),
    list("sinR_regulatory_K19449","sinR (Biofilm reg)",FALSE))
)
TD <- "phagcn3_high_apg/trees_v3"; HD <- "phagcn3_high_apg/hgt_tips_v3_strict"

for (cls in names(groups)) {
  plots <- lapply(groups[[cls]], function(it)
    make_circ(paste0(TD,"/",it[[1]],"/tree.treefile"),
              paste0(HD,"/",it[[1]],".txt"), it[[2]], it[[3]]))
  combined <- wrap_plots(plots, ncol=length(plots)) +
    plot_annotation(title=paste0(cls," AVGs — red = high-confidence HGT"),
                    theme=theme(plot.title=element_text(size=12,face="bold",hjust=0.5)))
  out <- paste0("hgtstrict_",tolower(cls),".png")
  ggsave(out, combined, width=5.5*length(plots), height=6.8, dpi=145)
  cat("saved",out,"\n")
}
cat("done\n")
