library(ggplot2)
library(reshape2)
library(viridis)
library(grid)
library(gtable)
library(RColorBrewer)
library(strucchange)
library(plyr)
library(ggrepel)
library(geosphere)
library("gridExtra")
library(measurements)
library(ggpubr)
library(tidyr)
library(scales)
library(grid)
library(measurements)
library(ggpubr)
library(plyr)
library(eulerr)
library(pvclust)
library("FactoMineR")
library("factoextra")
library(viridis)

#############
### FIGURE 1
#############

setwd("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/figure1/")
info <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/github/ressource/TEsClassificationLabeling.txt", h=F)
colnames(info) = c("idx", "sf", "sf_short_1", "sf_short_2", "sf_short_3")

#####
## FIGURE 1A: Ae. aegypti libraries composition
#####
mat <- read.table("all_libraries.classification.tab", h=F)
colnames(mat) = c("fam", "sf", "libname")
mat = merge(mat, info, by="sf")
ddply(mat, .(libname), plyr::summarize, count=length(fam))

mat = ddply(mat, .(libname, sf_short_2), plyr::summarize, count=length(fam))
mat$sf_short_2 = factor(mat$sf_short_2, levels = c("Unknown", "Other", "MITE", "TIR", "LINE", "LTR"))
mat$libname = factor(mat$libname, levels = c("denovo", "aag2", "Lau", "Tefam"))

cols <- c(LTR = "#1D78AE", LINE = "#2FAF8F", TIR = "#8B2F2F", MITE = "#C94C4C", Other = "#F2A3A3", Unknown = "#D6D6D6")

p1 = ggplot(mat, aes(x=libname, y=count, color=sf_short_2, fill=sf_short_2)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(legend.position="none") +
  ylab("# consensus")

p1

#ggsave("Aeaegypti.libComp.nbcons.pdf", width = 55, height = 55, units = "mm")

#####
## FIGURE 1B: Ae. aegypti libraries RM fraction
#####
classif <- read.table("all_libraries.classification.tab", h=F)
colnames(classif) = c("fam", "sf", "libname")

mat <- read.table("all_libraries.rm.tab", h=F)
colnames(mat) = c("fam", "count", "len", "libname")
mat = merge(mat, classif, by=c("fam", "libname"))
mat = merge(mat, info, by="sf")
ddply(mat, .(libname), plyr::summarize, len=sum(len)/1278732104)

mat = ddply(mat, .(libname, sf_short_2), plyr::summarize, len=sum(len))
total = ddply(mat, .(libname), plyr::summarize, len=sum(len))
total$len/1278732104
mat$sf_short_2 = factor(mat$sf_short_2, levels = c("Unknown", "Other", "MITE", "TIR", "LINE", "LTR"))
mat$libname = factor(mat$libname, levels = c("denovo", "aag2", "Lau", "Tefam"))
head(mat)

cols <- c(LTR = "#1D78AE", LINE = "#2FAF8F", TIR = "#8B2F2F", MITE = "#C94C4C", Other = "#F2A3A3", Unknown = "#D6D6D6")

p1 = ggplot(mat, aes(x=libname, y=len/1278732104*100, color=sf_short_2, fill=sf_short_2)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  ylim(c(0,100)) + theme(legend.position="none") +
  ylab("TE fraction (%)")
p1

# ggsave("Aeaegypti.libComp.TEfraction.pdf", width = 55, height = 55, units = "mm")

#####
## FIGURE 1D: LINE TE_1400_R1
#####

#### before curation
ncol <- 293 # TE_1400_R1 before
align <- read.table("TE_1400_R1.before.align.tab", h=F, fill = T, col.names=paste0('V', seq_len(ncol)))

# get order by start
order = paste(align[,1][order(align[,2], decreasing=TRUE)], align[,2][order(align[,2], decreasing=TRUE)], sep = "_")
order = unique(order)

align.df = data.frame()
for(i in 1:dim(align)[1]){
  n=length(align[i,!is.na(align[i,])])
  value = as.numeric(as.character(align[i,4:n]))
  pos = seq(align[i,2],align[i,2]+length(value)*10,10)
  length(pos)
  length(value)
  align[i,]
  
  df = data.frame(
    seq = paste(align[i,1], align[i,2], sep = "_"),
    pos = pos[1:length(value)],
    value = value
  )
  align.df = rbind(align.df, df)
}

align.df$seq = factor(align.df$seq, levels = order)
align.df$value = factor(align.df$value)

col = c("#ff6600", "#ffcc00", "#ccff00", "#66ff00", "#00ff00", "#00ff66", "#00ffcc", "#00ccff", "#0066ff", "#0000ff")
names(col) = 1:10

ggplot(data=align.df, aes(x=pos, y=seq, color=value, fill=value)) +
  geom_tile(height=0.25, width=0.25, size=0.25) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y = element_blank(),
    legend.position="none"
  )
# before 0.15 0.15 0.5
#ggsave("TE_1400_R1.before.pdf", width = 60, height = 50, units = "mm")

#### after curation
ncol <- 421 # TE_1400_R1 after
align <- read.table("TE_1400_R1.after.align.tab", h=F, fill = T, col.names=paste0('V', seq_len(ncol)))

# get order by start
order = paste(align[,1][order(align[,2], decreasing=TRUE)], align[,2][order(align[,2], decreasing=TRUE)], sep = "_")
order = unique(order)

align.df = data.frame()
for(i in 1:dim(align)[1]){
  n=length(align[i,!is.na(align[i,])])
  value = as.numeric(as.character(align[i,4:n]))
  pos = seq(align[i,2],align[i,2]+length(value)*10,10)
  length(pos)
  length(value)
  align[i,]
  
  df = data.frame(
    seq = paste(align[i,1], align[i,2], sep = "_"),
    pos = pos[1:length(value)],
    value = value
  )
  align.df = rbind(align.df, df)
}

align.df$seq = factor(align.df$seq, levels = order)
align.df$value = factor(align.df$value)

col = c("#ff6600", "#ffcc00", "#ccff00", "#66ff00", "#00ff00", "#00ff66", "#00ffcc", "#00ccff", "#0066ff", "#0000ff")
names(col) = 1:10

ggplot(data=align.df, aes(x=pos, y=seq, color=value, fill=value)) +
  geom_tile(height=0.25, width=0.25, size=0.25) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y = element_blank(),
    legend.position="none"
  )

# ggsave("TE_1400_R1.after.pdf", width = 60, height = 50, units = "mm")

#####
## FIGURE 1E: TE_0261_BelPao
#####

#### before curation
ncol <- 526 # TE_0261_BelPao.before
align <- read.table("TE_0261_BelPao.before.align.txt", h=F, fill = T, col.names=paste0('V', seq_len(ncol)))

# get order by start
order = paste(align[,1][order(align[,2], decreasing=TRUE)], align[,2][order(align[,2], decreasing=TRUE)], sep = "_")
order = unique(order)

align.df = data.frame()
for(i in 1:dim(align)[1]){
  n=length(align[i,!is.na(align[i,])])
  value = as.numeric(as.character(align[i,4:n]))
  pos = seq(align[i,2],align[i,2]+length(value)*10,10)
  length(pos)
  length(value)
  align[i,]
  
  df = data.frame(
    seq = paste(align[i,1], align[i,2], sep = "_"),
    pos = pos[1:length(value)],
    value = value
  )
  align.df = rbind(align.df, df)
}

align.df$seq = factor(align.df$seq, levels = order)
align.df$value = factor(align.df$value)

ggplot(data=align.df, aes(x=pos, y=seq, color=value, fill=value)) +
  geom_tile(height=0.25, width=0.25, size=0.15) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y = element_blank(),
    legend.position="none"
  )

ggsave("TE_0261_BelPao.before.pdf", width = 100, height = 100, units = "mm")

#### after curation
ncol <- 696 # TE_0261_BelPao.after
align <- read.table("TE_0261_BelPao.after.align.txt", h=F, fill = T, col.names=paste0('V', seq_len(ncol)))

# get order by start
order = paste(align[,1][order(align[,2], decreasing=TRUE)], align[,2][order(align[,2], decreasing=TRUE)], sep = "_")
order = unique(order)

align.df = data.frame()
for(i in 1:dim(align)[1]){
  n=length(align[i,!is.na(align[i,])])
  value = as.numeric(as.character(align[i,4:n]))
  pos = seq(align[i,2],align[i,2]+length(value)*10,10)
  length(pos)
  length(value)
  align[i,]
  
  df = data.frame(
    seq = paste(align[i,1], align[i,2], sep = "_"),
    pos = pos[1:length(value)],
    value = value
  )
  align.df = rbind(align.df, df)
}

align.df$seq = factor(align.df$seq, levels = order)
align.df$value = factor(align.df$value)

ggplot(data=align.df, aes(x=pos, y=seq, color=value, fill=value)) +
  geom_tile(height=0.25, width=0.25, size=0.15) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y = element_blank(),
    legend.position="none"
  )

#ggsave("TE_0261_BelPao.after.pdf", width = 100, height = 100, units = "mm")

#####
## Figure F: Assess completness of consensus
#####
mat_init <- read.table("aaegypti.ConsTElib.completness.txt", h=F)
colnames(mat_init) = c("consid", "sf", "status")
mat_init$version = "initial"
dim(mat_init)

mat_cur <- read.table("aaegypti.ConsTElib.currated.completness.txt", h=F)
colnames(mat_cur) = c("consid", "sf", "status")
mat_cur$version = "currated"

info <- read.table("/Volumes/@IVI/Josquin_Daron/TEsgenomics/paper/figures/ressource/TEsClassificationLabeling.txt", h=F)
colnames(info) = c("idx", "sf", "sf_short_1", "sf_short_2", "sf_short_3")

mat = rbind(mat_init[,c("consid", "sf", "status", "version")], mat_cur[,c("consid", "sf", "status", "version")])
head(mat)

mat = merge(mat, info, by="sf")
mat$sf_short_2[which(mat$sf_short_1=="Gypsy")] = "Gypsy"
mat$sf_short_2[which(mat$sf_short_1=="Copia")] = "Copia"
mat$sf_short_2[which(mat$sf_short_1=="BelPao")] = "BelPao"
head(mat)

# add remove tag 
mat = reshape2::dcast(mat, consid+sf+sf_short_2~version, value.var = "status")
mat$currated[is.na(mat$currated)] = "remove"
mat = mat[!is.na(mat$initial),] # remove TE that have change of superfamily names
mat = reshape2::melt(mat, measure.vars=c("currated", "initial"))
colnames(mat) = c("consid", "sf", "sf_short_2", "version", "status")
head(mat)

mat = ddply(mat, .(sf_short_2, status, version), plyr::summarize, count=length(consid))
mat$status = factor(mat$status, levels = c( "remove", "fragmented", "complete"))
dim(mat)

cols <- c(
  complete   = "#5B8DB8",  # steel blue
  fragmented = "#E8A838",  # amber
  remove     = "#BEBEBE"   # darker gray
)


ggplot(mat, aes(x=sf_short_2, y=count, color=NA, fill=status)) +
  geom_bar(stat="identity", position = "stack", width = 0.6) +
  facet_wrap(~version) +
  theme_classic() + 
  #  ggtitle("Library curration") +
  xlab("Superfamilies") +
  ylab("# of consensus") +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#    legend.position="none"
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)

#ggsave("Aeaegypti.libCompleteness.pdf", width = 240, height = 120, units = "mm")

# get stats
table(mat_cur$status)[1]/dim(mat_cur)[1]*100 # overall proportion of complete consensi in curated library

mat = rbind(mat_init[,c("consid", "sf", "status", "version")], mat_cur[,c("consid", "sf", "status", "version")])
mat = merge(mat, info, by="sf")
mat$sf_short_2[which(mat$sf_short_1=="Gypsy")] = "Gypsy"
mat$sf_short_2[which(mat$sf_short_1=="Copia")] = "Copia"
mat$sf_short_2[which(mat$sf_short_1=="BelPao")] = "BelPao"
mat = reshape2::dcast(mat, consid+sf+sf_short_2~version, value.var = "status")
mat$currated[is.na(mat$currated)] = "remove"
mat = mat[!is.na(mat$initial),] # remove TE that have change of superfamily names
length(which(mat$currated=="remove"))/dim(mat)[1]
