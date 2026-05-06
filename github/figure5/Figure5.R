library("FactoMineR")
library("factoextra")
library("gridExtra")
library(RColorBrewer)
library(eulerr)
library(geosphere)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggstance)
library(grid)
library(gtable)
library(measurements)
library(plyr)
library(pvclust)
library(reshape2)
library(scales)
library(strucchange)
library(tidyr)
library(viridis)
library(phytools)
library(phangorn)
library(purrr)
library(tidyr)
library(dplyr)
library(minpack.lm)  # for nlsLM (more robust than nls)
library(phytools)
library(phangorn)
library(ggtree)
library(ggtreeExtra)
library(data.table)
library(readr)
library(cowplot)
library(ape)
library(data.table)
library(utils)
library(scales)
library(beanplot)
library(GenomicRanges)
library(IRanges)
library(rstatix)
library(stringr)
library(ggforce)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggnewscale)
library(ggforce)
library(viridis)

#############
### FIGURE 5
#############

#####
## metadata
#####
info <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/TEsClassificationLabeling.txt", h=F)
colnames(info) = c("idx", "sf", "sf_short_1", "sf_short_2", "sf_short_3")

dmel <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/dmel.TEsClassification.tab", h=F)
colnames(dmel) = c("fam", "status", "sf")
dmel$sp = "dmel"

acoln1 <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/acoln3.TEsClassification.tab", h=F)
colnames(acoln1) = c("fam", "status", "sf")
acoln1$sp = "acoln3"

aaegl <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/aaegL5.TEsClassification.tab", h=F)
colnames(aaegl) = c("fam", "status", "sf")
aaegl = aaegl[which(aaegl$status!="remove"),]
aaegl$sp = "aaegL5"

classif = rbind(dmel[,c("fam", "sf", "status", "sp")], acoln1[,c("fam", "sf", "status", "sp")], aaegl[,c("fam", "sf", "status", "sp")])
classif = merge(classif, info, by="sf")

#####
## Figure 5 A: Example Dicer2 locus
#####

# shared sequence order / vertical positions
seq_ids <- c("AaegL5_3", "NC_064671.1", "2R")
y_map <- c("AaegL5_3" = 2, "NC_064671.1" = 1, "2R" = 0)

#--------------------------
# 1. Exons
#--------------------------
exon_cols <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

exon2 <- read.table(
  "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/figure5/exon.annotation.tab",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = ""
)
colnames(exon2) <- exon_cols

exon2 <- exon2 %>%
  dplyr::filter(type == "exon", seqid %in% seq_ids) %>%
  dplyr::mutate(
    Parent = stringr::str_match(attributes, "Parent=([^;]+)")[, 2],
    y = unname(y_map[seqid]),
    xmin = pmin(start, end),
    xmax = pmax(start, end),
    ymin = y - 0.08,
    ymax = y + 0.08
  ) %>%
  tidyr::separate_rows(Parent, sep = ",")

make_exon_arcs <- function(exon_df, arc_height = 0.12) {
  links <- exon_df %>%
    dplyr::arrange(seqid, Parent, xmin, xmax) %>%
    dplyr::group_by(seqid, Parent) %>%
    dplyr::mutate(
      next_xmin = dplyr::lead(xmin),
      this_xmax = xmax
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(next_xmin))
  
  do.call(rbind, lapply(seq_len(nrow(links)), function(i) {
    row <- links[i, ]
    x1 <- row$this_xmax
    x2 <- row$next_xmin
    xm <- (x1 + x2) / 2
    y0 <- row$y
    
    data.frame(
      arc_id = i,
      seqid = row$seqid,
      Parent = row$Parent,
      x = c(x1, xm, x2),
      y = c(y0, y0 + arc_height, y0)
    )
  }))
}

exon_arcs <- make_exon_arcs(exon2, arc_height = 0.12)

transcript_labels <- exon2 %>%
  dplyr::arrange(seqid, Parent, xmin) %>%
  dplyr::group_by(seqid, Parent) %>%
  dplyr::summarise(
    x = dplyr::first(xmin),
    y = dplyr::first(y) + 0.18,
    label = dplyr::first(Parent),
    .groups = "drop"
  )

#--------------------------
# 2. TEs
#--------------------------
te_cols <- c("seqid", "start", "end", "feature_id", "strand", "sf")

te <- read.table(
  "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/figures/figure5/TEs.annotation.bed",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)
colnames(te) <- te_cols

te <- merge(te, classif, by = "sf")

te_rects <- te %>%
  dplyr::filter(seqid %in% seq_ids) %>%
  dplyr::mutate(
    y = unname(y_map[seqid]),
    xmin = pmin(start, end),
    xmax = pmax(start, end),
    ymin = y - 0.10,
    ymax = y + 0.10
  )

#--------------------------
# 3. BLAST
#--------------------------
blast_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore")

b1 <- read.table(
  "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/figure5/aaegL5_vs_acoln3.tbltx",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)
b2 <- read.table(
  "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/figure5/acoln3_vs_dmel.tbltx",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)

colnames(b1) <- blast_cols
colnames(b2) <- blast_cols

prepare_blast <- function(df) {
  df %>%
    dplyr::mutate(
      qmin = pmin(qstart, qend),
      qmax = pmax(qstart, qend),
      smin = pmin(sstart, send),
      smax = pmax(sstart, send),
      strand = ifelse((qend - qstart) * (send - sstart) > 0, "forward", "reverse")
    )
}

b1 <- prepare_blast(b1)
b2 <- prepare_blast(b2)

make_poly <- function(df, y_top, y_bottom, id_offset = 0, h = 0.14) {
  do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    x <- df[i, ]
    data.frame(
      hit = id_offset + i,
      x = c(x$qmin, x$qmax, x$smax, x$smin),
      y = c(y_top - h, y_top - h, y_bottom + h, y_bottom + h),
      pident = x$pident,
      bitscore = x$bitscore,
      strand = x$strand
    )
  }))
}

poly1 <- make_poly(b1, y_top = 2, y_bottom = 1, id_offset = 0, h = 0.14)
poly2 <- make_poly(b2, y_top = 1, y_bottom = 0, id_offset = 100000, h = 0.14)
poly_all <- dplyr::bind_rows(poly1, poly2)

seq_df <- data.frame(
  seqname = c(unique(b1$qseqid), unique(b1$sseqid), unique(b2$sseqid)),
  y = c(2, 1, 0),
  xmin = c(min(b1$qmin), min(b1$smin), min(b2$smin)),
  xmax = c(max(b1$qmax), max(b1$smax), max(b2$smax))
)

#--------------------------
# 4. Colors
#--------------------------
te_cols_fill <- c(
  Gypsy   = "#1D78AE",
  BelPao  = "#37a2e0",
  Copia   = "#8FD0F2",
  LINE    = "#2FAF8F",
  TIR     = "#8B2F2F",
  MITE    = "#C94C4C",
  Other   = "#F2A3A3",
  Unknown = "#D6D6D6"
)

#--------------------------
# 5. Plot
#--------------------------
ggplot() +
  geom_polygon(
    data = poly_all,
    aes(x = x, y = y, group = hit, fill = pident),
    alpha = 0.55, color = NA
  ) +
  scale_fill_viridis_c(name = "% identity") +
  ggnewscale::new_scale_fill() +
  geom_segment(
    data = seq_df,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    linewidth = 1, color = "black"
  ) +
  geom_rect(
    data = te_rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = sf_short_3),
    color = NA
  ) +
  geom_rect(
    data = exon2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white", color = "darkblue", linewidth = 0.3
  ) +
  ggforce::geom_bezier(
    data = exon_arcs,
    aes(x = x, y = y, group = arc_id),
    color = "darkblue", linewidth = 0.4
  ) +
  # Uncomment if you want transcript labels
  # geom_text(
  #   data = transcript_labels,
  #   aes(x = x, y = y, label = label),
  #   hjust = 0, vjust = 0,
  #   size = 3, color = "darkblue"
  # ) +
  geom_text(
    data = seq_df,
    aes(x = xmin, y = y + 0.15, label = seqname),
    hjust = 0, size = 4
  ) +
  scale_fill_manual(values = te_cols_fill, name = "TE class") +
  theme_bw() +
  labs(x = "Genomic position", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ggsave("figure4A.TEatDicer2Locus.pdf", units = "mm", width = 200, height = 70)

#####
## Figure 5 C: Distance genes to nearest TE
#####

# classification info
info <- fread("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/TEsClassificationLabeling.txt", header = FALSE)
setnames(info, c("idx", "sf", "sf_short_1", "sf_short_2", "sf_short_3"))

dmel <- fread("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/dmel.TEsClassification.tab", header = FALSE)
setnames(dmel, c("fam", "status", "sf"))
dmel[, sp := "dmel"]

acoln1 <- fread("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/acoln3.TEsClassification.tab", header = FALSE)
setnames(acoln1, c("fam", "status", "sf"))
acoln1[, sp := "acoln3"]

aaegl <- fread("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/aaegL5.TEsClassification.tab", header = FALSE)
setnames(aaegl, c("fam", "status", "sf"))
aaegl <- aaegl[status != "remove"]
aaegl[, sp := "aaegL5"]

classif_dt <- rbindlist(list(
  dmel[, .(fam, sf, status, sp)],
  acoln1[, .(fam, sf, status, sp)],
  aaegl[, .(fam, sf, status, sp)]
))

classif_dt <- merge(classif_dt, info, by = "sf")
classif_dt <- unique(classif_dt, by = c("sf", "sp"))

# distance files
files <- c(
  aaegL5 = "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/figure5/aaegL5.GenesclosedTEs.txt.gz",
  acoln3 = "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/figure5/acoln3.GenesclosedTEs.txt.gz",
  dmel   = "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/figure5/dmel.GenesclosedTEs.txt.gz"
)

dist_list <- lapply(names(files), function(sp) {
  dt <- fread(files[[sp]], sep = "\t", header = FALSE)
  setnames(dt, c("id", "strand", "sf", "fam", "distance"))
  dt[, sp := sp]
  dt
})

dist <- rbindlist(dist_list)
dist <- dist[distance > 0]

# join using sf + sp
dist[classif_dt, sf_short_3 := i.sf_short_3, on = c("sf", "sp")]
dist$sf_short_3[which(dist$sf_short_3%in%c("Gypsy", "Copia", "BelPao"))] = "LTR"
dist$sf_short_3 = factor(dist$sf_short_3, levels = c("LTR", "LINE", "MITE", "TIR", "Other", "Unknown"))

stat_df <- dist %>%
  group_by(sf_short_3) %>%
  kruskal_test(distance ~ sp) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(
    label = case_when(
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01  ~ "**",
      p.adj <= 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

y_pos <- dist %>%
  dplyr::group_by(sf_short_3) %>%
  dplyr::summarise(y = max(distance, na.rm = TRUE) * 1.5)

stat_df <- left_join(stat_df, y_pos, by = "sf_short_3")

ggplot(dist, aes(x = sf_short_3, y = distance, color = sp)) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(
    data = stat_df,
    aes(x = sf_short_3, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000, 10000, 100000),
    labels = label_number()
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    #    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor.y = element_blank()
  ) +
  labs(x = "TE class", y = "Distance to nearest gene")

ddply(dist, .(sp), plyr::summarize, avg=mean(distance))

# ggsave("figure5B.distanceG2T.pdf", units = "mm", width = 90, height = 60)

#####
## Figure 5 D: TE content in introns
#####

species_cols <- c(
  aaegL5 = "#6A4C93",  # muted purple
  acoln3 = "#E69F00",  # warm orange
  dmel   = "#4D4D4D"   # dark gray
)

# genome size
genomeSize = data.frame(
  sp = c("dmel", "acoln3", "aaegL5"),
  bp_total_genome = c(143726002,249560815,1278732104)
)

info <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/TEsClassificationLabeling.txt", h=F)
colnames(info) = c("idx", "sf", "sf_short_1", "sf_short_2", "sf_short_3")

dmel <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/dmel.TEsClassification.tab", h=F)
colnames(dmel) = c("fam", "status", "sf")
dmel$sp = "dmel"

acoln1 <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/acoln3.TEsClassification.tab", h=F)
colnames(acoln1) = c("fam", "status", "sf")
acoln1$sp = "acoln3"

aaegl <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/aaegL5.TEsClassification.tab", h=F)
colnames(aaegl) = c("fam", "status", "sf")
aaegl = aaegl[which(aaegl$status!="remove"),]
aaegl$sp = "aaegL5"

classif = rbind(dmel[,c("fam", "sf", "status", "sp")], acoln1[,c("fam", "sf", "status", "sp")], aaegl[,c("fam", "sf", "status", "sp")])
classif = merge(classif, info, by="sf")

# load intron bed files
files <- c("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/figures/figure5/aaegL5_intron.v2.bed",
           "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/figures/figure5/acoln3_intron.v2.bed",
           "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/figures/figure5/dmel_intron.v2.bed")

list <- lapply(files, function(fn) {
  read.table(fn, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
})

intron <- do.call(rbind, list)
colnames(intron) = c("chrom", "start", "end", "id", "score", "strand", "sp")

# load TEs annotation bed files
files <- c("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/aaegL5.TEs.bed",
           "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/acoln3.TEs.bed",
           "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/dmel.TEs.bed")

list <- lapply(files, function(fn) {
  read.table(fn, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
})

mat <- do.call(rbind, list)
colnames(mat) = c("chrom", "start", "end", "teid", "strand", "sf", "fam")

mat = merge(mat, classif, by=c("sf" ,"fam"))

# find overlap
gr_te <- makeGRangesFromDataFrame(
  mat,
  seqnames.field="chrom", start.field="start", end.field="end",
  keep.extra.columns=TRUE
)

gr_intron <- makeGRangesFromDataFrame(
  intron,
  seqnames.field="chrom", start.field="start", end.field="end",
  keep.extra.columns=TRUE
)

# overlaps TE x intron
hits <- findOverlaps(gr_te, gr_intron)

mat$in_intron <- "F"
mat$in_intron[unique(queryHits(hits))] <- "T"
mat$sf_short_3[which(mat$sf_short_3%in%c("Gypsy", "Copia", "BelPao"))] = "LTR"
mat.sum = ddply(mat, .(sp, sf_short_3, in_intron), plyr::summarize, bp=sum(end-start+1))
mat.total = ddply(intron, .(sp), plyr::summarize, bp_total_intron=sum(end-start+1))
mat.sum = merge(mat.sum, mat.total, by="sp")
mat.sum = mat.sum[which(mat.sum$in_intron=="T"),]
mat.sum$perc = mat.sum$bp/mat.sum$bp_total_intron*100

mat.sum$sf_short_3 = factor(mat.sum$sf_short_3, levels = c("LTR", "LINE", "MITE", "TIR", "Other", "Unknown"))
mat.sum = rbind(mat.sum, data.frame(sp="dmel", sf_short_3="MITE", in_intron="T", bp=NA, bp_total_intron=NA, perc=0))

ggplot(mat.sum, aes(x = sf_short_3, y = perc, color = sp, fill = sp)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# ggsave("figure5C.TEcontentInIntrons.pdf", units = "mm", width = 90, height = 60)

