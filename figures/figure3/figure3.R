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
### FIGURE 3
#############

#####
## metadata
#####
setwd("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/figure3/")

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
## Figure 3A: K2P distance to the consensus
#####
step = 1
start <- seq(0, 100, step)
end <- seq(step, 100+step, step)
intv <- data.frame(cbind(start, end))

zz=gzfile('dmel.consDistance.txt.gz','rt')
dmel=read.csv(zz,header=F,sep = "\t")
colnames(dmel) = c("chrom", "start", "end", "fam", "startCons", "endCons", "id")

dmel$intv = NA
for(i in 1:dim(intv)[1]){
  dmel$intv[which(dmel$id>=intv$start[i] & dmel$id<intv$end[i])] = intv$start[i]
}

zz=gzfile('acoln3.consDistance.txt.gz','rt')
acoln3=read.csv(zz,header=F,sep = "\t")
colnames(acoln3) = c("chrom", "start", "end", "fam", "startCons", "endCons", "id")
head(acoln3)

acoln3$intv = NA
for(i in 1:dim(intv)[1]){
  acoln3$intv[which(acoln3$id>=intv$start[i] & acoln3$id<intv$end[i])] = intv$start[i]
}

zz=gzfile('aaegL5.consDistance.txt.gz','rt')
aaegL5=read.csv(zz,header=F,sep = "\t")
colnames(aaegL5) = c("chrom", "start", "end", "fam", "startCons", "endCons", "id")
head(aaegL5)

sum(aaegL5$end-aaegL5$start+1, na.rm = TRUE)

aaegL5$intv = NA
for(i in 1:dim(intv)[1]){
  aaegL5$intv[which(aaegL5$id>=intv$start[i] & aaegL5$id<intv$end[i])] = intv$start[i]
}

# merge and calculation
mat = rbind(dmel, acoln3, aaegL5)
mat = ddply(mat, .(fam, intv), plyr::summarize, len=sum(end-start+1))
mat = merge(mat, classif, by="fam")
mat = ddply(mat, .(sf_short_3, intv, sp), plyr::summarize, len=sum(len))
mat = merge(mat, genomeSize, "sp")

mat$sf_short_3[which(mat$sf_short_3%in%c("Gypsy", "BelPao", "Copia"))] = "LTR"
mat$sf_short_3 = factor(mat$sf_short_3, levels = rev(c("LTR", "LINE", "TIR", "MITE", "Other", "Unknown")))

# create dummy 
dummy = data.frame()
for(sp in unique(mat$sp)){
  for(sf in unique(mat$sf_short_3)){
    tmp = intv
    tmp$sp=sp
    tmp$sf_short_3=sf
    tmp$len = 0
    tmp$size = genomeSize$size[which(genomeSize$sp==sp)]
    dummy = rbind(dummy, tmp)
  }
}

colnames(dummy)[which(colnames(dummy)=="end")] = "intv"
head(dummy)

# merge dummy and mat
mat = rbind(mat, dummy[,c("sp", "sf_short_3", "intv", "len", "size")])
mat = ddply(mat, .(sp, sf_short_3, intv, size), plyr::summarize, len=sum(len))

cols <- c(
  LTR            = "#1D78AE",
  LINE           = "#2FAF8F",
  TIR            = "#8B2F2F",
  MITE           = "#C94C4C",
  Other          = "#F2A3A3",
  Unknown        = "#D6D6D6"
)

ggplot(mat, aes(x=intv, y=len/size*100, fill=sf_short_3)) +
  geom_bar(stat="identity", position = "stack") +
  facet_wrap(~sp, ncol = 1) +
  theme_classic() +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  xlim(c(-1,60)) +
  theme(
    legend.position="none",
    strip.background = element_blank()
  ) +
  xlab("Kimura divergence from consensus (CpG-adjused)") +
  ylab("Genome fraction (%)")

# ggsave("TEs_dynamics.stack_sf.K2P.pdf", units = "mm", width = 100, height = 125)

sum(mat$len[which(mat$sp=="dmel" & mat$intv==0)])/143726002
sum(mat$len[which(mat$sp=="acoln3" & mat$intv==0)])/249560815
sum(mat$len[which(mat$sp=="aaegL5" & mat$intv==0)])/1278732104

#####
## Figure 3B: Estimation decay rate using distance to consensus
#####
mat$dist = mat$len/mat$size*100

# get rate
eDecay <- function(t, ampl, tau) ampl * exp(-t / tau)

fit_one <- function(dat) {
  dat <- dat %>%
    transmute(t = intv, y = dist) %>%
    filter(is.finite(t), is.finite(y)) %>%
    arrange(t)
  
  if (nrow(dat) < 3 || all(dat$y <= 0)) return(NULL)
  
  A0 <- max(dat$y, na.rm = TRUE)
  t_half <- dat$t[which(dat$y < A0 / 2)][1]
  tau0 <- if (!is.na(t_half) && is.finite(t_half) && t_half > 0) t_half else median(dat$t[dat$t > 0], na.rm = TRUE)
  
  fit <- tryCatch(
    nlsLM(y ~ eDecay(t, ampl, tau),
          data = dat,
          start = list(ampl = A0, tau = tau0),
          lower = c(ampl = 0, tau = 1e-12),
          control = nls.lm.control(maxiter = 200)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  co <- coef(fit)
  
  # SE from covariance matrix (approximate, Wald-type)
  V <- tryCatch(vcov(fit), error = function(e) NULL)
  tau_se  <- if (!is.null(V) && "tau" %in% rownames(V)) sqrt(V["tau","tau"]) else NA_real_
  ampl_se <- if (!is.null(V) && "ampl" %in% rownames(V)) sqrt(V["ampl","ampl"]) else NA_real_
  
  tibble(
    ampl = unname(co["ampl"]),
    tau  = unname(co["tau"]),
    ampl_se = ampl_se,
    tau_se  = tau_se
  )
}

fits <- mat %>%
  group_by(sp, sf_short_3) %>%
  group_modify(~{
    out <- fit_one(.x)
    if (is.null(out)) tibble(ampl = NA_real_, tau = NA_real_, ampl_se = NA_real_, tau_se = NA_real_) else out
  }) %>%
  ungroup()

# CI + transform to rate scale
df_tau <- fits %>%
  filter(is.finite(tau)) %>%
  mutate(
    tau_lwr = pmax(tau - 1.96 * tau_se, 1e-12),
    tau_upr = tau + 1.96 * tau_se,
    rate    = 1 / tau,
    rate_lwr = 1 / tau_upr,
    rate_upr = 1 / tau_lwr
  ) %>%
  select(sp, sf_short_3, tau, tau_se, tau_lwr, tau_upr, rate, rate_lwr, rate_upr)

#df_tau$sf <- factor(df_tau$sf, levels = c("Gypsy","BelPao","Copia","LINE","TIR"))
df_tau = rbind(df_tau, data.frame(sp="dmel", sf_short_3="MITE", tau=NA, tau_se=NA, tau_lwr=NA, tau_upr=NA, rate=NA, rate_lwr=NA, rate_upr=NA))
df_tau = df_tau[which(df_tau$sf_short_3%in%c("LTR","LINE","TIR","MITE")),]
df_tau$sf_short_3 <- factor(df_tau$sf_short_3, levels = c("LTR","LINE","TIR","MITE"))
df_tau$sp <- factor(df_tau$sp)

pd <- position_dodge(width = 0.8)

ggplot(df_tau, aes(x = sp, y = rate, fill = sf_short_3)) +
  geom_col(aes(group = sf_short_3), position = pd, width = 0.7) +
  geom_errorbar(aes(ymin = rate_lwr, ymax = rate_upr, group = sf_short_3),
                position = pd, width = 0.2, linewidth = 0.4) +
  theme_classic() +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = expression("TE decay rate (1/" * tau * ")")) +
  theme(
    legend.position="none",
    strip.background = element_blank()
  )

# ggsave("TEs_dynamics.slope.disttocons.pdf", units = "mm", width = 100, height = 40)

#####
## Figure 3C: TE dynamics using Terminal branch length distance
#####
dmel_term=read.table("dmel.termBranchLen.tab",header=F,sep = "\t")
colnames(dmel_term) = c("chrom", "start", "end", "fam", "len", "sp", "sf_short_3", "clstr")

acol_term=read.table("acoln3.termBranchLen.tab",header=F,sep = "\t")
colnames(acol_term) = c("chrom", "start", "end", "fam", "len", "sp", "sf_short_3", "clstr")

aaegl5_term=read.table("aaegL5.termBranchLen.tab",header=F,sep = "\t")
colnames(aaegl5_term) = c("chrom", "start", "end", "fam", "len", "sp", "sf_short_3", "clstr")

mat = rbind(dmel_term, acol_term, aaegl5_term)
mat$sf_short_3[which(mat$sf_short_3%in%c("Sola", "Tc1-Mariner"))] = "TIR"
mat = merge(mat, classif, id=c("fam", "sp", "sf_short_3"))
mat$sf_short_3[which(mat$sf_short_3%in%c("BelPao", "Copia", "Gypsy"))] = "LTR"

# get distribution
df = data.frame()
for(sp in unique(mat$sp)){
  for(sf in unique(mat$sf_short_3)){
    sub = mat[which(mat$sp==sp & mat$sf_short_3==sf),]
    x <- sub$len[sub$len >= 0 & sub$len <= 0.05]
    h = hist(x, breaks = seq(0,0.05,0.001), plot = FALSE)
    tmp = data.frame(breaks=seq(0,0.05,0.001)[1-(length(h$counts)-1)], count=h$counts, sp=sp, sf=sf)
    df = rbind(df, tmp)
  }
}

#df$sf = factor(df$sf, levels = rev(c("Gypsy", "BelPao", "Copia", "LINE", "TIR")))
df$sf = factor(df$sf, levels = rev(c("LTR", "LINE", "TIR")))

cols <- c(Gypsy = "#1D78AE", BelPao = "#37a2e0", Copia = "#8FD0F2", LINE = "#2FAF8F", TIR = "#8B2F2F")
cols <- c(LTR = "#1D78AE", LINE = "#2FAF8F", TIR = "#8B2F2F")

df = merge(df, genomeSize, by="sp")

ggplot() +
  geom_bar(data=df, aes(x=breaks, y=count/size*1e6, color=sf, fill=sf), stat = "identity", position = "stack") +
  #  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = c("1", "10", "100", "1000", "10000", "100000")) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  facet_wrap(~sp, ncol=1) +
  theme_classic() +
  theme(
    legend.position="none",
    strip.background = element_blank()
  ) +
  labs(y="Number of TE copies", x="Terminal bramch length distance (sub/site)")

# ggsave("TEs_dynamics.stack.termBranch.pdf", units = "mm", width = 100, height = 125)

sum(df$count[which(df$sp=="dmel" & df$breaks==0)])/143726002*1e6
sum(df$count[which(df$sp=="acoln3" & df$breaks==0)])/249560815*1e6
sum(df$count[which(df$sp=="aaegL5" & df$breaks==0)])/1278732104*1e6

#####
## Figure 3D: Estimation decay rate using terminal branch length distance
#####

# get rate
eDecay <- function(t, ampl, tau) ampl * exp(-t / tau)

fit_one <- function(dat) {
  dat <- dat %>%
    transmute(t = breaks, y = count) %>%
    filter(is.finite(t), is.finite(y)) %>%
    arrange(t)
  
  if (nrow(dat) < 3 || all(dat$y <= 0)) return(NULL)
  
  A0 <- max(dat$y, na.rm = TRUE)
  t_half <- dat$t[which(dat$y < A0 / 2)][1]
  tau0 <- if (!is.na(t_half) && is.finite(t_half) && t_half > 0) t_half else median(dat$t[dat$t > 0], na.rm = TRUE)
  
  fit <- tryCatch(
    nlsLM(y ~ eDecay(t, ampl, tau),
          data = dat,
          start = list(ampl = A0, tau = tau0),
          lower = c(ampl = 0, tau = 1e-12),
          control = nls.lm.control(maxiter = 200)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  co <- coef(fit)
  
  # SE from covariance matrix (approximate, Wald-type)
  V <- tryCatch(vcov(fit), error = function(e) NULL)
  tau_se  <- if (!is.null(V) && "tau" %in% rownames(V)) sqrt(V["tau","tau"]) else NA_real_
  ampl_se <- if (!is.null(V) && "ampl" %in% rownames(V)) sqrt(V["ampl","ampl"]) else NA_real_
  
  tibble(
    ampl = unname(co["ampl"]),
    tau  = unname(co["tau"]),
    ampl_se = ampl_se,
    tau_se  = tau_se
  )
}

fits <- df %>%
  group_by(sp, sf) %>%
  group_modify(~{
    out <- fit_one(.x)
    if (is.null(out)) tibble(ampl = NA_real_, tau = NA_real_, ampl_se = NA_real_, tau_se = NA_real_) else out
  }) %>%
  ungroup()

# CI + transform to rate scale
df_tau <- fits %>%
  filter(is.finite(tau)) %>%
  mutate(
    tau_lwr = pmax(tau - 1.96 * tau_se, 1e-12),
    tau_upr = tau + 1.96 * tau_se,
    rate    = 1 / tau,
    rate_lwr = 1 / tau_upr,
    rate_upr = 1 / tau_lwr
  ) %>%
  select(sp, sf, tau, tau_se, tau_lwr, tau_upr, rate, rate_lwr, rate_upr)

#df_tau$sf <- factor(df_tau$sf, levels = c("Gypsy","BelPao","Copia","LINE","TIR"))
df_tau$sf <- factor(df_tau$sf, levels = c("LTR","LINE","TIR"))
df_tau$sp <- factor(df_tau$sp)

pd <- position_dodge(width = 0.8)

ggplot(df_tau, aes(x = sp, y = rate, fill = sf)) +
  geom_col(aes(group = sf), position = pd, width = 0.7) +
  geom_errorbar(aes(ymin = rate_lwr, ymax = rate_upr, group = sf),
                position = pd, width = 0.2, linewidth = 0.4) +
  theme_classic() +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = expression("TE decay rate (1/" * tau * ")")) +
  theme(
    legend.position="none",
    strip.background = element_blank()
  )

#ggsave("TEs_dynamics.slope.termBranch.pdf", units = "mm", width = 100, height = 40)

