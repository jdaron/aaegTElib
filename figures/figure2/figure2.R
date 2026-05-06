library("FactoMineR")
library("factoextra")
library("gridExtra")
library(RColorBrewer)
library(data.table)
library(dplyr)
library(eulerr)
library(geosphere)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggstance)
library(ggtree)
library(ggtreeExtra)
library(grid)
library(gtable)
library(measurements)
library(minpack.lm)  # for nlsLM (more robust than nls)
library(phangorn)
library(phytools)
library(plyr)
library(purrr)
library(pvclust)
library(readr)
library(reshape2)
library(scales)
library(strucchange)
library(tidyr)
library(viridis)

#############
### FIGURE 2
#############

#####
## metadata
#####
setwd("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/figure2/")

info <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/TEsClassificationLabeling.txt", h=F)
colnames(info) = c("idx", "sf", "sf_short_1", "sf_short_2", "sf_short_3")

dmel <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/dmel.TEsClassification.tab", h=F)
colnames(dmel) = c("fam", "status", "sf")
dmel$sp = "dmel"

acoln3 <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/acoln3.TEsClassification.tab", h=F)
colnames(acoln3) = c("fam", "status", "sf")
acoln3$sp = "acoln3"

aaegl <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/aaegL5.TEsClassification.tab", h=F)
colnames(aaegl) = c("fam", "status", "sf")
aaegl$sp = "aaegL5"

classif = rbind(dmel[,c("fam", "sf", "status", "sp")], acoln3[,c("fam", "sf", "status", "sp")], aaegl[,c("fam", "sf", "status", "sp")])
classif = merge(classif, info, by="sf")
head(classif)

# genome size
genomeSize = data.frame(
  sp = c("dmel", "acoln3", "aaegL5"),
  size = c(143726002,249560815,1278732104)
)

#####
## Figure 2A: TE fraction donuts
#####
files <- c("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/aaegL5.TEs.bed",
           "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/acoln3.TEs.bed",
           "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/dmel.TEs.bed")

list <- lapply(files, function(fn) {
  read.table(gzfile(fn, "rt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
})

mat <- do.call(rbind, list)
colnames(mat) <- c("chrom", "start", "end", "teid", "strand", "sf", "fam")
mat = merge(mat, classif[,c("sf", "fam", "sp", "sf_short_3")], by=c("sf", "fam"))
genomefraction = ddply(mat, .(sf_short_3, sp), plyr::summarize, len=sum(end-start))
exon = data.frame(sf_short_3=c("exon", "exon", "exon"), sp=c("dmel", "acoln3", "aaegL5"), len=c(36066829, 44334954, 37682663))
genomefraction = rbind(genomefraction, exon)
genomefraction = merge(genomefraction, genomeSize, by="sp")
head(genomefraction)

tmp = ddply(genomefraction, .(sp, size), plyr::summarize, whole=sum(len))
tmp$len = tmp$size-tmp$whole
tmp$sf_short_3 = "low_complexity"

genomefraction = rbind(genomefraction, tmp[,c("sp", "sf_short_3", "len", "size")])
genomefraction$perc = genomefraction$len/genomefraction$size*100

genomefraction$sf_short_3 = factor(genomefraction$sf_short_3, levels = c("low_complexity", "exon", "Unknown", "Other", "MITE", "TIR", "LINE", "Copia", "BelPao", "Gypsy"))

cols <- c(
  low_complexity = "#E3B44D",
  exon           = "#9B8FD1",
  LTR            = "#1D78AE",
  Gypsy          = "#1D78AE",
  BelPao         = "#37a2e0",
  Copia          = "#8FD0F2",
  LINE           = "#2FAF8F",
  TIR            = "#8B2F2F",
  MITE           = "#C94C4C",
  Other          = "#F2A3A3",
  Unknown        = "#D6D6D6"
)

ggplot(genomefraction, aes(x = 2, y = perc, fill = sf_short_3)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +                      # controls hole size
  facet_wrap(~ sp, nrow = 1) +          # three facets for three species
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = NULL, fill = "sf_short_3") +
  theme_void() +
  scale_fill_manual(values = cols) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

# ggsave("TEsFractioninGenome.pdf", units = "mm", width = 100, height = 70)
ddply(genomefraction[which(genomefraction$sf_short_3%in%c("Gypsy", "Copia", "BelPao", "LINE", "TIR", "MITE", "Other", "Unknown")),], .(sp), plyr::summarize, perc=sum(perc))

#####
## Figure 2B: Boxplot TE fragmentation
#####
len <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/aaegL5_acoln3_dmel.consensusLength.txt", h=F)
colnames(len) = c("sp", "fam", "len")
classif_len = merge(classif, len, by=c("sp", "fam"))

files <- c("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/aaegL5.complete_fragment.tab.gz",
           "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/acoln3.complete_fragment.tab.gz",
           "/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/dmel.complete_fragment.tab.gz")

list <- lapply(files, function(fn) {
  read.table(gzfile(fn, "rt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
})

mat <- do.call(rbind, list)
colnames(mat) <- c("chrom", "start", "end", "sf", "fam", "startc", "endc")

mat = merge(mat, classif_len, by=c("sf", "fam"))
mat$perc = (mat$endc-mat$startc+1)/mat$len*100

mat = mat[which(mat$status=="complete"),]
mat = mat[which(mat$sf_short_3%in%c("BelPao", "Copia", "Gypsy", "LINE", "TIR", "MITE")),]
mat$sf_short_3[which(mat$sf_short_3%in%c("BelPao", "Copia", "Gypsy"))] = "LTR"
mat$sf_short_3 = factor(mat$sf_short_3, c("LTR", "LINE", "TIR", "MITE"))

ggplot(mat, aes(x = sf_short_3, y = perc, color = sf_short_3)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~sp, scales = "free_x", ncol=1) +
  scale_color_manual(values = cols) +
  theme_bw() +
  theme(legend.position="none") +
  labs(x="", y="TE copy length (% of consensus length)")  +
  ylim(c(0,100))

# ggsave("TEfragments.boxplot.pdf", units = "mm", width = 60, height = 130)

# stats
library(dplyr)
library(emmeans)   # post‑hoc

head(mat[,c("sp", "sf_short_3", "perc")])
mat <- mat %>%
  mutate(sp = factor(sp),
         sf_short_3 = factor(sf_short_3))

# include interaction
fit <- aov(perc ~ sp * sf_short_3, data = mat)
summary(fit)

# Post‑hoc
emmeans(fit, ~ sp | sf_short_3) %>% pairs(adjust = "tukey")      # compare species within each TE class
emmeans(fit, ~ sf_short_3 | sp) %>% pairs(adjust = "tukey")      # compare TE classes within each species

#####
## Diversity barplot
#####
mat <- read.table("all3sp.nbCopiesPerFam.denovoClassif.txt", h=F)
colnames(mat) = c("sp", "sf", "fam", "count")

mat$sf[which(mat$sf=="Sola")] = "TIR"
mat$sf[which(mat$sf=="Tc1-Mariner")] = "TIR"
mat$sf[which(mat$sf%in%c("BelPao", "Copia", "Gypsy"))] = "LTR"

mat <- mat %>%
  dplyr::group_by(sp, sf) %>%
  dplyr::arrange(desc(count), fam, .by_group = TRUE) %>%  # sort by count within sp×sf
  dplyr::mutate(fam_idx = row_number()) %>%               # rank index within sp×sf
  dplyr::ungroup()

mat = mat[which(mat$count>=3),]

mat.d = ddply(mat, .(sp, sf), plyr::summarize, count=length(count))
#cols <- c(Gypsy = "#1D78AE", BelPao = "#37a2e0", Copia = "#8FD0F2", LINE = "#2FAF8F", TIR = "#8B2F2F")
cols <- c(LTR = "#1D78AE", LINE = "#2FAF8F", TIR = "#8B2F2F")
mat$sf = factor(mat$sf, levels=c("LTR", "LINE", "TIR"))

ggplot(mat %>% filter(count > 0, fam_idx > 0),
       aes(fam_idx, count, color = sf, fill = sf)) +
  geom_col() +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  facet_grid(sp ~ sf) +
  scale_y_log10(
    breaks = breaks_log(n = 6, base = 10),
    labels = label_number()
  ) +
  labs(x = "fam_idx", y = "Count") +
  theme_bw() +
  theme(
    legend.position = "none"
    #    strip.background = element_blank()
  ) +
  labs(x="TE families (ranked by copy number)", y="Number of TE copies per family")

# ggsave("diversityTEs.pdf", units = "mm", width = 80, height = 110)

