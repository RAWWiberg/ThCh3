##
# D.pseudoobscura project evolution linse
# PoPoolation2/glm data analysis
# Last Modified: January 2017
##
#clean environment
rm(list = ls(all = TRUE))
#
# Load libraries and define functions
# #####
library(plyr)
library(ggplot2)
library(scales)
#install.packages("grid")
library(gridExtra)
library(grid)
#install.packages("reshape")
library(reshape)
library(dplyr)
library(qvalue)
se <- function(x) {sqrt(var(x, na.rm = TRUE))/sqrt(length(x))}
source("~/Desktop/Data/RData/RScripts/ggplot_theme.R")

# ####
setwd("~/Desktop/Data/RData/pseudo_evol/")
#
# Assess Genome Coverage
# ####
# Load the data
R1M1_cov_hist<-read.table("R1M1_genome_covhist.txt",
           header = FALSE)
# Subset to coverage < 400 and > 0
colnames(R1M1_cov_hist)<-c("level","coverage","count","total","prop")
R1M1_cov_hist <- R1M1_cov_hist[R1M1_cov_hist$coverage < 400 &
                                 R1M1_cov_hist$coverage > 0,]
R1M1_cov_hist$data<-rep("R1M1",nrow(R1M1_cov_hist))
R1M1_cov_hist$cumprop <- cumsum(R1M1_cov_hist$prop)
max(R1M1_cov_hist$coverage[R1M1_cov_hist$cumprop < 0.5])

R1P1_cov_hist<-read.table("R1P1_genome_covhist.txt",
  header = FALSE)
# Subset to coverage < 400 and > 0
colnames(R1P1_cov_hist)<-c("level","coverage","count","total","prop")
R1P1_cov_hist <- R1P1_cov_hist[R1P1_cov_hist$coverage < 400 &
                                 R1P1_cov_hist$coverage > 0,]
R1P1_cov_hist$data<-rep("R1P1",nrow(R1P1_cov_hist))
R1P1_cov_hist$cumprop <- cumsum(R1P1_cov_hist$prop)
max(R1P1_cov_hist$coverage[R1P1_cov_hist$cumprop < 0.5])

R2M2_cov_hist<-read.table("R2M2_genome_covhist.txt",
  header = FALSE)
# Subset to coverage < 400 and > 0
colnames(R2M2_cov_hist)<-c("level","coverage","count","total","prop")
R2M2_cov_hist <- R2M2_cov_hist[R2M2_cov_hist$coverage < 400 &
                                 R2M2_cov_hist$coverage > 0,]
R2M2_cov_hist$data<-rep("R2M2",nrow(R2M2_cov_hist))
R2M2_cov_hist$cumprop <- cumsum(R2M2_cov_hist$prop)
max(R2M2_cov_hist$coverage[R2M2_cov_hist$cumprop < 0.5])

R2P1_cov_hist<-read.table("R2P1_genome_covhist.txt",
  header = FALSE)
# Subset to coverage < 400 and > 0
colnames(R2P1_cov_hist)<-c("level","coverage","count","total","prop")
R2P1_cov_hist <- R2P1_cov_hist[R2P1_cov_hist$coverage < 400 &
                                 R2P1_cov_hist$coverage > 0,]
R2P1_cov_hist$data<-rep("R2P1",nrow(R2P1_cov_hist))
R2P1_cov_hist$cumprop <- cumsum(R2P1_cov_hist$prop)
max(R2P1_cov_hist$coverage[R2P1_cov_hist$cumprop < 0.5])

R3M1_cov_hist<-read.table("R3M1_genome_covhist.txt",
  header = FALSE)
# Subset to coverage < 400 and > 0
colnames(R3M1_cov_hist)<-c("level","coverage","count","total","prop")
R3M1_cov_hist <- R3M1_cov_hist[R3M1_cov_hist$coverage < 400 &
                                 R3M1_cov_hist$coverage > 0,]
R3M1_cov_hist$data<-rep("R3M1",nrow(R3M1_cov_hist))
R3M1_cov_hist$cumprop <- cumsum(R3M1_cov_hist$prop)
max(R3M1_cov_hist$coverage[R3M1_cov_hist$cumprop < 0.5])

R3P2_cov_hist<-read.table("R3P2_genome_covhist.txt",
  header = FALSE)
# Subset to coverage < 400 and > 0
colnames(R3P2_cov_hist)<-c("level","coverage","count","total","prop")
R3P2_cov_hist <- R3P2_cov_hist[R3P2_cov_hist$coverage < 400 &
                                 R3P2_cov_hist$coverage > 0,]
R3P2_cov_hist$data<-rep("R3P2",nrow(R3P2_cov_hist))
R3P2_cov_hist$cumprop <- cumsum(R3P2_cov_hist$prop)
max(R3P2_cov_hist$coverage[R3P2_cov_hist$cumprop < 0.5])

R4M1_cov_hist<-read.table("R4M1_genome_covhist.txt",
  header = FALSE)
# Subset to coverage < 400 and > 0
colnames(R4M1_cov_hist)<-c("level","coverage","count","total","prop")
R4M1_cov_hist <- R4M1_cov_hist[R4M1_cov_hist$coverage < 400 &
                                 R4M1_cov_hist$coverage > 0,]
R4M1_cov_hist$data<-rep("R4M1",nrow(R4M1_cov_hist))
R4M1_cov_hist$cumprop <- cumsum(R4M1_cov_hist$prop)
max(R4M1_cov_hist$coverage[R4M1_cov_hist$cumprop < 0.5])

R4P1_cov_hist<-read.table("R4P1_genome_covhist.txt",
  header = FALSE)
# Subset to coverage < 400 and > 0
colnames(R4P1_cov_hist)<-c("level","coverage","count","total","prop")
R4P1_cov_hist <- R4P1_cov_hist[R4P1_cov_hist$coverage < 400 &
                                 R4P1_cov_hist$coverage > 0,]
R4P1_cov_hist$data<-rep("R4P1",nrow(R4P1_cov_hist))
R4P1_cov_hist$cumprop <- cumsum(R4P1_cov_hist$prop)
max(R4P1_cov_hist$coverage[R4P1_cov_hist$cumprop < 0.5])

#combine the data
all_cov <- rbind(R1M1_cov_hist,
                 R1P1_cov_hist,
                R2M2_cov_hist,
                R2P1_cov_hist,
                R3M1_cov_hist,
                R3P2_cov_hist,
                R4M1_cov_hist,
                R4P1_cov_hist)
# What is the mean coverage
head(all_cov)

str(all_cov)
head(all_cov)
all_cov$data <- as.factor(all_cov$data)
unique(all_cov$data)
# Plot the histogram
cov_hist<-ggplot()+
  geom_bar(data=all_cov[all_cov$coverage > 0,],
           aes(x=coverage,y=prop),
           stat="identity")+
#  scale_x_continuous(breaks=seq(0,400,10))+
  geom_vline(xintercept = c(17,48),colour="red",linetype="dashed")+
  ylab("Proportion of Sites")+
  xlab("Coverage")+
  facet_grid(data~.)

cov_hist+my.theme+theme(axis.text.x = element_text(size = 12,
                                                   angle=45,
                                                   vjust = 1,
                                                   hjust = 1),
                        axis.text.y = element_text(size = 12),
                        strip.text.y=element_text(size=12,angle=0))

# Coverage threshholds of 30 and 400 are kind of arbitrary but probably
# sensible.

# Get percentiles for the aggregated distribution
aggr_cov<-aggregate(all_cov$count,list(coverage=all_cov$coverage),sum)
aggr_cov$prop <- aggr_cov$x/sum(as.numeric(aggr_cov$x))
aggr_cov$data<-rep("aggregate",nrow(aggr_cov))

colnames(aggr_cov)<-c("coverage","count","prop","data")

aggr_cov$prop2<-aggr_cov$count/sum(as.numeric(aggr_cov$count))
aggr_cov$cumprop <- cumsum(aggr_cov$prop2)

# 10th percentile
max(aggr_cov$coverage[aggr_cov$cumprop < 0.1])

# 90th percentile
max(aggr_cov$coverage[aggr_cov$cumprop < 0.9])

head(aggr_cov[c(1,2,3,4)])
head(all_cov[c(2,3,5,6)])
all_cov <- rbind(all_cov[c(2,3,5,6)],
                 aggr_cov[c(1,2,3,4)])
tail(all_cov)

vlin<-data.frame(xint1=17,xint2=49,data2="Aggregate")
unique(all_cov$data2)
all_cov$data2<-factor(all_cov$data,
                     labels=c("M Line\nReplicate 1",
                              "E Line\nReplicate 1",
                              "M Line\nReplicate 2",
                              "E Line\nReplicate 2",
                              "M Line\nReplicate 3",
                              "E Line\nReplicate 3",
                              "M Line\nReplicate 4",
                              "E Line\nReplicate 4",
                              "Aggregate"))
cov_hist<-ggplot()+
  geom_bar(data=all_cov[all_cov$coverage > 0 & all_cov$coverage < 100,],
           aes(x=coverage,y=prop),
           stat="identity")+
  #  scale_x_continuous(breaks=seq(0,400,10))+
  geom_vline(data=vlin,aes(xintercept = xint1),
             colour="red",linetype="dashed")+
  geom_vline(data=vlin,aes(xintercept = xint2),colour="red",linetype="dashed")+
  scale_y_continuous(breaks=c(0,0.02,0.04))+
  ylab("Proportion of Sites")+
  xlab("Coverage")+
  facet_grid(data2~.)

cov_hist+my.theme+theme(axis.text.x = element_text(size = 10,
                                                   angle=45,
                                                   vjust = 1,
                                                   hjust = 1),
                        axis.title = element_text(size=10),
                        axis.text.y = element_text(size = 10),
                        strip.text.y=element_text(size=12,angle=0))

# Which percentiles correspond to coverage 30 and 400
# ####

#
# Standard CMH-test results.
# ####
# Load data
# min cov = 17, max cov = 49
cmhdat <- read.table("pseudo_evol_minc16_cov17-49.cmh",
                     sep = "", header = FALSE)
colnames(cmhdat) <- c("Chr","Pos","ref",
                      "R1M","R1P","R2M","R2P",
                      "R3M","R3P","R4M","R4P",
                      "cmh_p")
cmhdat$SNP <- paste("SNP_",cmhdat$Chr,"_",cmhdat$Pos, sep="")
# Get top 1% p-value
cmh_top1 <- quantile(cmhdat$cmh_p, 0.01, na.rm = TRUE)
cmh_top1
# Get Bonferroni threshhold
cmh_bonf<-0.05/nrow(cmhdat)

# Subset top 1%
cmhdat_top<-cmhdat[which(cmhdat$cmh_p < cmh_top1),]
cmhdat_top$SNPtype <- rep("Top 1%",nrow(cmhdat_top))

# Background SNPs (everything not in topSNPs)
cmhdat_bcgr <- cmhdat[which(!(cmhdat$SNP %in% cmhdat_top$SNP)),]
cmhdat_bcgr$SNPtype <- rep("Background",nrow(cmhdat_bcgr))

# Combine
cmhdat <- rbind(cmhdat_top, cmhdat_bcgr)

#subset to known chromosomes (for plotting)
cmhdat_chr <- cmhdat[which(substr(cmhdat$Chr,1,1) != "U"),]
cmhdat_chr$Chr <- as.character(cmhdat_chr$Chr)

nrow(cmhdat_chr[cmhdat_chr$SNPtype == "Top 1%",])
nrow(cmhdat_chr[cmhdat_chr$cmh_p < cmh_bonf,])
nrow(cmhdat_chr[cmhdat_chr$SNPtype == "Background",])

#order chromosome parts 
cmhdat_chr$Chr2 <- factor(gsub("_group",".",cmhdat_chr$Chr))
cmhdat_chr <- sort_df(cmhdat_chr, vars = c("Chr2","Pos"))
cmhdat_chr$Chr <- factor(gsub("_.*","",cmhdat_chr$Chr)) 

#make nice vectors
cmhdat_chr$absPos <- seq(1:nrow(cmhdat_chr))
cmhdat_chr$Chr <- factor(cmhdat_chr$Chr)
chrNum=length(unique(cmhdat_chr$Chr))
c <- 1
midvec <- vector(length = chrNum)
for (i in unique(cmhdat_chr$Chr)){
  ndx <- which(cmhdat_chr$Chr == i)
  SubPos <- cmhdat_chr$absPos[ndx]
  midvec[c] <-((max(SubPos) - min(SubPos))/2) + min(SubPos)
  c <- c + 1
}

cmh_manhplot <- ggplot() +
  geom_point(data = cmhdat_chr[cmhdat_chr$SNPtype != "Top 1%",], 
             aes(Pos/1000000,
                 -log10(cmh_p), colour = Chr),alpha = 1/5) +
  geom_point(data = cmhdat_chr[cmhdat_chr$SNPtype == "Top 1%",], 
             aes(Pos/1000000,
                 -log10(cmh_p)),colour = "red",alpha = 1/5) +
  geom_hline(yintercept = -log10(cmh_bonf),
             linetype = "dashed",colour="red")+
  #  geom_hline(yintercept = -log10(glm_top5),
  #             linetype = "dashed",colour="orange",size = 1)+
  #  geom_hline(yintercept = -log10(glm_top10),
  #             linetype = "dashed",colour="blue",size = 1)+
  scale_colour_manual(limits = as.character(unique(cmhdat_chr$Chr)),
                      values = c(rep(c(
                        "grey50","black"),2),"grey50"), 
                      guide = FALSE) +
  #  scale_x_continuous(labels=unique(cmhdat_chr$Chr), 
  #                     breaks = midvec) +
  scale_y_continuous(limits=c(0,75),breaks = c(0,30,60)) +
  xlab("Position on Chromosome (Mb)") +
  ylab("-log10(CMH p-value)")+
  facet_grid(Chr2~.)

cmh_manhplot + my.theme  +
  theme(axis.text.x  = element_text(size = 12, 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        axis.text.y  = element_text(size = 12, 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        strip.text.y = element_text(size = 10,angle = 0))

# ####



# QUASIBINOMIAL GLM results.
# ####
# Load Data
#sc<-"_neff"
sc<-""
# min cov = 17, max cov = 49, scale = raw counts
glmdat <- read.table(paste("pseudo_evol_qbglm-unp",
                           sc,
                           "_minc16_cov17-49.rout",sep=""),
                     sep = "", header = FALSE)
colnames(glmdat) <- c("Chr","Pos","ref",
                      "R1M","R1P","R2M","R2P",
                      "R3M","R3P","R4M","R4P",
                      "t_p")
head(glmdat)

# t = Treatment (E vs M)[Line]
glmdat$SNP <- paste("SNP_",glmdat$Chr,"_",glmdat$Pos, sep="")
glmdat$Chr <- gsub("^", "chrom", glmdat$Chr)
glmdat$Chr2 <- factor(gsub("_group",".",glmdat$Chr))

head(glmdat)
nrow(glmdat)
str(glmdat)


#####
##Without subsetting to known chromosomes for getting top 1%
#####

#find top 1% of SNPs
#distribution of p-values
glm_t_p_hist <- ggplot() +
  geom_histogram(data = glmdat, aes(t_p)) +
  xlab("p-value") +
  ylab("")
glm_t_p_hist<-glm_t_p_hist + my.theme+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_blank())
# Q-Q plot
o<-sort(glmdat$t_p,decreasing=FALSE)
pdat<-data.frame(o=o,
                 e=(1:length(o))/length(o))
head(pdat)
nrow(pdat)
q_q_plot<-ggplot()+
  geom_point(data=pdat,aes(-log10(e),-log10(o)),size=0.2)+
  geom_line(data=pdat,aes(-log10(e),-log10(e)),
            linetype="dashed",colour="red")+
  xlab("Expected -log10(p-values)")+
  ylab("Observed -log10(p-values)")+
  my.theme+
  theme(axis.text = element_text(size=10),
        axis.title.y = element_text(size=14))

# Make Q-Q plot with p-value distribution inset.
xleft   = 0.05
xright  = 0.45
ybottom = 0.55
ytop    = 0.95 

# Calculate position in plot1 coordinates
# Extract x and y values from q_q_plot
l1 = ggplot_build(q_q_plot)
x1 = l1$layout$panel_ranges[[1]]$x.range[1]
x2 = l1$layout$panel_ranges[[1]]$x.range[2]
y1 = l1$layout$panel_ranges[[1]]$y.range[1]
y2 = l1$layout$panel_ranges[[1]]$y.range[2]
xdif = x2-x1
ydif = y2-y1
xmin  = x1 + (xleft*xdif)
xmax  = x1 + (xright*xdif)
ymin  = y1 + (ybottom*ydif)
ymax  = y1 + (ytop*ydif) 

# Get plot2 and make grob
g2 = ggplotGrob(glm_t_p_hist)
plot3 = q_q_plot + annotation_custom(grob = g2, 
                                     xmin=xmin, xmax=xmax, 
                                     ymin=ymin, ymax=ymax)
plot(plot3)

glm_t_p_hist <- ggplot() +
  geom_histogram(data = glmdat, aes(t_p)) +
  xlab("T p-value") +
  geom_vline(xintercept = 0.05)+
  ylab("Count")
glm_t_p_hist + my.theme

#####

# Find Top SNPS
# ####

#TopSNPs
# Get SNPs that have q-values < 0.05
# Save t_p as list
write.table(glmdat$t_p,paste("glmdat",sc,"_t_p.list",sep=""),
            row.names = FALSE,col.names = FALSE,quote=FALSE)
# Read the p-values
pvals<-scan(paste("glmdat",sc,"_t_p.list",sep=""))
# Convert to q-values with default settings
qobj<-qvalue(pvals)
qplot(qobj)
# Attach q-values to data
glmdat$t_q<-qobj$qvalues
# Which SNPs have t_q < 0.05
glmdat_top_q<-glmdat[
  glmdat$t_q < 0.05,]
glmdat_top_q$SNPtype <- rep("q-value",nrow(glmdat_top_q))
nrow(glmdat_top_q)
#Background SNPs (everything not in topSNPs)
glmdat_bcgr_q <- glmdat[which(!(glmdat$SNP %in% glmdat_top_q$SNP)),]
glmdat_bcgr_q$SNPtype <- rep("Background",nrow(glmdat_bcgr_q))


# Overlap in the top SNPs
head(glmdat_top_q)

#Combine
glmdat_q <- rbind(glmdat_top_q, glmdat_bcgr_q)
nrow(glmdat_q)
# Save the top q-value SNPs
write.table(glmdat_top_q,
            paste("pseudo_evol_qbglm-unp",
                  sc,
                  "_minc16_cov17-49_top_qval_snps.rout",sep=""),
            row.names = FALSE,col.names = TRUE,quote=FALSE)

#
## Manhattan plot
#####
head(glmdat_chr)
#subset to known chromosomes (for plotting)
glmdat_chr <- glmdat_q[grep("Unknown",glmdat_q$Chr,invert=TRUE),]
glmdat_chr$Chr <- as.character(glmdat_chr$Chr)
nrow(glmdat_chr[glmdat_chr$SNPtype == "q-value",])

#order chromosome parts 
glmdat_chr$Chr2 <- factor(gsub("chrom","",glmdat_chr$Chr2))
glmdat_chr <- sort_df(glmdat_chr, vars = c("Chr2","Pos"))
glmdat_chr$Chr <- factor(gsub("_.*","",glmdat_chr$Chr)) 
nrow(glmdat_chr)
#make nice vectors
glmdat_chr$absPos <- seq(1:nrow(glmdat_chr))
glmdat_chr$Chr <- factor(glmdat_chr$Chr)
chrNum=length(unique(glmdat_chr$Chr))
c <- 1
midvec <- vector(length = chrNum)
for (i in unique(glmdat_chr$Chr)){
  ndx <- which(glmdat_chr$Chr == i)
  SubPos <- glmdat_chr$absPos[ndx]
  midvec[c] <-((max(SubPos) - min(SubPos))/2) + min(SubPos)
  c <- c + 1
}
# Draw manhattan plot of -log10(T p-values)
# All snps
head(glmdat_chr)
unique(glmdat_chr$Chr)
unique(glmdat_chr$Chr2)
manhplot <- ggplot() +
  geom_point(data = glmdat_chr[glmdat_chr$SNPtype != "q-value",], 
             aes(Pos/1000000,
                 -log10(t_p)),colour = "grey70",alpha = 1/5) +
  geom_point(data = glmdat_chr[glmdat_chr$SNPtype == "q-value",], 
             aes(Pos/1000000,
                 -log10(t_p)),colour = "red",alpha = 1/5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",colour="red")+
  scale_x_continuous(breaks=c(0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30),
                     labels=c("0","","5","","10","","15","","20","","25","","30"))+
  xlab("Position on Chromosome (Mb)") +
  ylab("-log10(T p-value)")+
  facet_grid(Chr2~.)
manhplot + my.theme  +
#  scale_y_continuous(limits=c(0,6),breaks = c(0,2.5,5)) +
  theme(axis.text.x  = element_text(size = 10, 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        axis.text.y  = element_text(size = 10, 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        axis.title = element_text(size = 10),
        strip.text.y = element_text(size = 10,angle = 0))


# Make region dataset that highlights the "chimney-peaks"
chim_reg_dat<-data.frame(Chr2=c("2",
                                "3",
                                "4.3","4.4",
                                "XL.1e","XL.1e",
                                "XR.6","XR.6",
                                "XR.8","XR.8"),
                         reg=c("1",
                               "2",
                               "3","4",
                               "5","6",
                               "7","8",
                               "9","10"),
                     start=c(5200000,
                             18750000,
                             3700000,2500000,
                             1250000,9500000,
                             5600000,8500000,
                             1300000,2800000),
                     end=c(6700000,
                           19800000,
                           4500000,3700000,
                           2500000,11500000,
                           6500000,10000000,
                           2000000,4000000))
manhplot + 
  geom_rect(data=chim_reg_dat,aes(xmin=start/1000000,xmax=end/1000000,
                                  ymin=4.8,ymax=5.4))+
  geom_text(data=chim_reg_dat,aes((end+200000)/1000000,y=4.7,label=reg),size=4)+
  my.theme  +
  scale_y_continuous(limits=c(0,6),breaks = c(0,5)) +
  theme(axis.text.x  = element_text(size = 12, 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        axis.text.y  = element_text(size = 10, 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        strip.text.y = element_text(size = 10,angle = 0))

# Breakpoint data from Wallace et al., 2011
breaks<-data.frame(Chr2=c("3","3","3","3","3",
                          "3","3","3","3",
                          "3","3"),
                   x=c(2292730,6432192,6472192,8900335,9140888,
                       9832232,10830454,14259585,15426271,
                       17444472,17705213),
                   name=c("pSTPP","pHYSC","pSCTL","pSTAR","pHYST",
                          "dSTPP","dSCTL","pSCCH","dSCCH",
                          "dHYSC","dHYST"))
# Draw just 3rd chromosome with breakpoints
chr3manhplot <- ggplot() +
  geom_point(data = glmdat_chr[glmdat_chr$SNPtype != "q-value" & 
                                 glmdat_chr$Chr2 == "3",], 
             aes(Pos/1000000,
                 -log10(t_q)), colour = "grey70",alpha = 1/5) +
  geom_point(data = glmdat_chr[glmdat_chr$SNPtype == "q-value" & 
                                 glmdat_chr$Chr2 == "3",], 
             aes(Pos/1000000,
                 -log10(t_q)),colour = "red",alpha = 1/5) +
  geom_vline(data=breaks,aes(xintercept = x/1000000),
             size=0.5,
             colour=c("red","purple","lightblue","darkblue",
                      "green","red","lightblue","darkorange",
                      "darkorange","purple","green"))+
  scale_y_continuous(limits=c(0,5),breaks = c(0,2.5,5)) +
  scale_x_continuous(breaks=c(0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30),
                     labels=c("0","","5","","10","","15","","20","","25","","30"))+
  xlab("Position on Chromosome (Mb)") +
  ylab("-log10(T q-value)")+
  facet_grid(Chr2~.)

chr3manhplot + my.theme  +
  theme(axis.text.x  = element_text(size = 12, 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        axis.text.y  = element_text(size = 12, 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        axis.title = element_text(size = 10),
        strip.text.y = element_text(size = 10,angle = 0))

# Distribution of SNPs:
# Chi squared test of nr of top SNPs on main chroms
unique(glmdat_chr$Chr)
# Expected distribution should be related to the proportional 
# length of each chromosome
chrom_dat_chr<-chrom_dat[grep("Unknown",chrom_dat$Chr,invert = TRUE),]
chrom_dat_chr$Chr2<-gsub("_group.*","",chrom_dat_chr$Chr)
chrom_dat_chr[,c(9,6)]
chr_l<-tapply(chrom_dat_chr$length,INDEX=list(chrom_dat_chr$Chr2),sum)
chr_l/sum(chr_l)

# All SNPS
tapply(glmdat_chr$Chr, INDEX=list(glmdat_chr$Chr),length)

chisq.test(tapply(glmdat_chr$Chr, INDEX=list(glmdat_chr$Chr),length),
           p=chr_l/sum(chr_l))
str(chisq.test(tapply(glmdat_chr$Chr, INDEX=list(glmdat_chr$Chr),length),
               p=chr_l/sum(chr_l)))
# Just "Background" SNPs
chisq.test(tapply(glmdat_chr$Chr[glmdat_chr$SNPtype=="Background"],
                  INDEX=list(glmdat_chr$Chr[glmdat_chr$SNPtype=="Background"]),
                  length),
           p=chr_l/sum(chr_l))
str(chisq.test(tapply(glmdat_chr$Chr[glmdat_chr$SNPtype=="Background"],
                  INDEX=list(glmdat_chr$Chr[glmdat_chr$SNPtype=="Background"]),
                  length),
           p=chr_l/sum(chr_l)))

# Top SNPs
tapply(glmdat_chr$Chr[glmdat_chr$SNPtype=="q-value"],
       INDEX=list(glmdat_chr$Chr[glmdat_chr$SNPtype=="q-value"]),
       length)
chisq.test(tapply(glmdat_chr$Chr[glmdat_chr$SNPtype=="q-value"],
                  INDEX=list(glmdat_chr$Chr[glmdat_chr$SNPtype=="q-value"]),
                  length),
           p=chr_l/sum(chr_l))
str(chisq.test(tapply(glmdat_chr$Chr[glmdat_chr$SNPtype=="q-value"],
                  INDEX=list(glmdat_chr$Chr[glmdat_chr$SNPtype=="q-value"]),
                  length),
           p=chr_l/sum(chr_l)))

#
# Which genes are between 18,750,000 and 19,750,000?
# ####

# Load annotation
dpse_ann<-read.table("~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/Dpse_genes_dmelnames.gtf",header = FALSE,sep="\t")

head(dpse_ann)
levels(dpse_ann$V1)
chim_reg_dat$region<-paste("region",seq(1,nrow(chim_reg_dat)),sep="")

chim_reg_ann<-data.frame()
for (region in chim_reg_dat$region){
    start<-chim_reg_dat$start[chim_reg_dat$region == region]
    end<-chim_reg_dat$end[chim_reg_dat$region == region]
    chr<-paste("chrom",
               chim_reg_dat$Chr[chim_reg_dat$region == region],
               sep="")
    chr <- gsub("\\.","_group",chr)
    
    dpse_ann_sub<-dpse_ann[
      dpse_ann$V1 == chr & dpse_ann$V4 >= start & dpse_ann$V5 <= end,
      c(1,4,5,9)]
    dpse_ann_sub$V9<-gsub("gene_id FBGN","FBgn",dpse_ann_sub$V9)
    dpse_ann_sub$V9<-gsub(";","",dpse_ann_sub$V9)
    dpse_ann_sub$V9<-gsub("_.","",dpse_ann_sub$V9)
    chim_reg_ann<-rbind(chim_reg_ann,dpse_ann_sub)
}
colnames(chim_reg_ann)<-c("Chr","start","end","gene")

head(chim_reg_ann)

# How many exons?
nrow(chim_reg_ann)
# Which genes
length(unique(chim_reg_ann$gene))

# Write the genes as a list
write.table(unique(chim_reg_ann$gene),"~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/dpse_glm_bonf_peaks_regions.list",
            quote=FALSE,row.names=FALSE,col.names=FALSE)

# Write the genes as a table
write.table(chim_reg_ann,"~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/dpse_glm_bonf_peaks_regions.tab",
            quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

# Load genes table
chim_reg_ann<-read.table("~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/dpse_glm_bonf_peaks_regions.tab",
            header=FALSE,sep="\t")
colnames(chim_reg_ann)<-c("Chr","start","end","FBid")

# Load gene ID conversion table
chim_reg_ann_conv<-read.table("~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/dpse_glm_bonf_peaks_regions_FBconv.list",
                         header=FALSE,sep="\t")
colnames(chim_reg_ann_conv)<-c("geneidsubmit","FBgeneid",
                          "FBgeneidconv","genename")
head(chim_reg_ann_conv)

# Add the gene name to the table
chim_reg_ann$name<-vector(length=nrow(chim_reg_ann))
for(id in chim_reg_ann$FBid){
  genename<-as.character(chim_reg_ann_conv$genename[
    chim_reg_ann_conv$geneidsubmit == id])
  chim_reg_ann$name[chim_reg_ann$FBid == id]<-genename  
}
head(chim_reg_ann)

# Write the genes as a table
write.table(chim_reg_ann,"~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/TableS3_dpse_glm_bonf_peaks_regions.tab",
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

# ####


#
# Write SNPs to files
# ####
# Load chromosome info
chrom_dat<-read.table("~/Desktop/Data/pseudo_project/reference/dpse_r31_FB2013_02/pseudo_chrom_data.tab",
                      header = TRUE, sep = ",")

#All SNPs
nrow(glmdat)
nrow(glmdat_bonf)
nrow(glmdat_q)

str(glmdat$Pos)

# Will use glmdat_bonf just for convenience
glmdat_1 <- glmdat_bonf
gtfdat <- cbind(glmdat_1$Chr3,
                rep("PoPoolation2",nrow(glmdat_1)),
                rep("SNP",nrow(glmdat_1)),
                glmdat_1$Pos,
                glmdat_1$Pos,
                rep(".",nrow(glmdat_1)),
                rep(".",nrow(glmdat_1)),
                rep(".",nrow(glmdat_1)),
                paste(
                  paste('gene_id "SNP_',
                        seq(1,nrow(glmdat_1),1),'"',sep=""),
                  paste('transcript_id "SNP_',
                        seq(1,nrow(glmdat_1),1),'"',sep=""),
                  paste('T-p_value "',glmdat_1$t_p,'"',sep=""),
                  sep = "; "))

gtfdat <- as.data.frame(gtfdat)
head(gtfdat)
nrow(gtfdat)

#SNP +- 30 bp
gtfdat30bp <- cbind(as.character(glmdat_1$Chr),
                    rep("PoPoolation2",nrow(glmdat_1)),
                    rep("CDS",nrow(glmdat_1)),
                    (glmdat_1$Pos-30),
                    (glmdat_1$Pos+30),
                    rep(".",nrow(glmdat_1)),
                    rep(".",nrow(glmdat_1)),
                    rep(".",nrow(glmdat_1)),
                    paste(
                      paste('gene_id "SNP_',
                            seq(1,nrow(glmdat_1),1),'"',sep=""),
                      paste('transcript_id "SNP_',
                            seq(1,nrow(glmdat_1),1),'"',sep=""),
                      paste('T-p_value "',
                            glmdat_1$t_p,'"',sep=""),
                      sep = "; "))

head(gtfdat30bp[,c(4,5)])
head(gtfdat[,c(4,5)])

gtfdat30bp<-as.data.frame(gtfdat30bp)
gtfdat30bp$V4 <- as.integer(as.character(gtfdat30bp$V4))
gtfdat30bp$V5 <- as.integer(as.character(gtfdat30bp$V5))
str(gtfdat30bp)
str(chrom_dat)
nrow(gtfdat30bp)
head(gtfdat30bp)
summary(chrom_dat)

# If any start positions < 0 set to 0
# If any end positions > max length set to max length
# This will stop regions going off the ends of chromosomes.
for (i in seq(1,nrow(gtfdat30bp))){
  if (gtfdat30bp$V4[i] < 0){
    gtfdat30bp$V4[i] <- 0
  }
  if (gtfdat30bp$V5[i] > chrom_dat$length[
    chrom_dat$Chr == as.character(gtfdat30bp$V1[i])]){
    gtfdat30bp$V5[i] <- chrom_dat$length[
      chrom_dat$Chr == as.character(gtfdat30bp$V1[i])]
  }
}
gtfdat30bp$V4 <- as.integer(as.character(gtfdat30bp$V4))
gtfdat30bp$V5 <- as.integer(as.character(gtfdat30bp$V5))

 
# Write top SNPs to a gtf file
head(glmdat_bonf)
head(glmdat_q)

top_bonfSNPs <- glmdat_bonf[glmdat_bonf$SNPtype == "Bonferroni",]
top_qSNPs <- glmdat_q[glmdat_q$SNPtype == "q-value",]
nrow(top_bonfSNPs)
nrow(top_qSNPs)

top_bonfSNPsgtfdat <- cbind(top_bonfSNPs$Chr3,
                       rep("PoPoolation2",nrow(top_bonfSNPs)),
                       rep("SNP",nrow(top_bonfSNPs)),
                       top_bonfSNPs$Pos,
                       top_bonfSNPs$Pos,
                       rep(".",nrow(top_bonfSNPs)),
                       rep(".",nrow(top_bonfSNPs)),
                       rep(".",nrow(top_bonfSNPs)),
                       paste(
                         paste('gene_id ','"',
                               top_bonfSNPs$SNP,'"',sep=""),
                         paste('transcript_id ','"',
                               top_bonfSNPs$SNP,'"',sep=""),
                         paste('T-p_value "',
                               top_bonfSNPs$t_p,'"',sep=""),
                         sep = "; "))
nrow(top_bonfSNPsgtfdat)

top_qSNPsgtfdat <- cbind(top_qSNPs$Chr3,
                       rep("PoPoolation2",nrow(top_qSNPs)),
                       rep("SNP",nrow(top_qSNPs)),
                       top_qSNPs$Pos,
                       top_qSNPs$Pos,
                       rep(".",nrow(top_qSNPs)),
                       rep(".",nrow(top_qSNPs)),
                       rep(".",nrow(top_qSNPs)),
                       paste(
                         paste('gene_id ','"',
                               top_qSNPs$SNP,'"',sep=""),
                         paste('transcript_id ','"',
                               top_qSNPs$SNP,'"',sep=""),
                         paste('T-p_value "',
                               top_qSNPs$t_p,'"',sep=""),
                         sep = "; "))
nrow(top_qSNPsgtfdat)

top_bonfSNPsgtfdat30bp <- cbind(as.character(top_bonfSNPs$Chr),
                           rep("PoPoolation2",nrow(top_bonfSNPs)),
                       rep("CDS",nrow(top_bonfSNPs)),
                       (top_bonfSNPs$Pos-30),
                       (top_bonfSNPs$Pos+30),
                       rep(".",nrow(top_bonfSNPs)),
                       rep(".",nrow(top_bonfSNPs)),
                       rep(".",nrow(top_bonfSNPs)),
                       paste(
                         paste('gene_id "SNP_',
                               seq(1,nrow(top_bonfSNPs),1),'"',sep=""),
                         paste('transcript_id "SNP_',
                               seq(1,nrow(top_bonfSNPs),1),'"',sep=""),
                         paste('T-p_value "',
                               top_bonfSNPs$t_p,'"',sep=""),
                         sep = "; "))


top_bonfSNPsgtfdat30bp<-as.data.frame(top_bonfSNPsgtfdat30bp)
top_bonfSNPsgtfdat30bp$V4 <- as.integer(as.character(
  top_bonfSNPsgtfdat30bp$V4))
top_bonfSNPsgtfdat30bp$V5 <- as.integer(as.character(
  top_bonfSNPsgtfdat30bp$V5))
# if any start positions < 0 set to 0
for (i in seq(1,nrow(top_bonfSNPsgtfdat30bp))){
  if (top_bonfSNPsgtfdat30bp$V4[i] < 0){
    top_bonfSNPsgtfdat30bp$V4[i] <- 0
  }
  if (top_bonfSNPsgtfdat30bp$V5[i] > chrom_dat$length[
    chrom_dat$Chr == as.character(top_bonfSNPsgtfdat30bp$V1[i])]){
    top_bonfSNPsgtfdat30bp$V5[i] <- chrom_dat$length[
      chrom_dat$Chr == as.character(top_bonfSNPsgtfdat30bp$V1[i])]
  }
}
nrow(top_bonfSNPsgtfdat30bp)
top_bonfSNPsgtfdat30bp$V4 <- as.integer(as.character(
  top_bonfSNPsgtfdat30bp$V4))
top_bonfSNPsgtfdat30bp$V5 <- as.integer(as.character(
  top_bonfSNPsgtfdat30bp$V5))

top_qSNPsgtfdat30bp <- cbind(as.character(top_qSNPs$Chr),
                           rep("PoPoolation2",nrow(top_qSNPs)),
                           rep("CDS",nrow(top_qSNPs)),
                           (top_qSNPs$Pos-30),
                           (top_qSNPs$Pos+30),
                           rep(".",nrow(top_qSNPs)),
                           rep(".",nrow(top_qSNPs)),
                           rep(".",nrow(top_qSNPs)),
                           paste(
                             paste('gene_id "SNP_',
                                   seq(1,nrow(top_qSNPs),1),'"',sep=""),
                             paste('transcript_id "SNP_',
                                   seq(1,nrow(top_qSNPs),1),'"',sep=""),
                             paste('T-p_value "',
                                   top_qSNPs$t_p,'"',sep=""),
                             sep = "; "))


top_qSNPsgtfdat30bp<-as.data.frame(top_qSNPsgtfdat30bp)
top_qSNPsgtfdat30bp$V4 <- as.integer(as.character(top_qSNPsgtfdat30bp$V4))
top_qSNPsgtfdat30bp$V5 <- as.integer(as.character(top_qSNPsgtfdat30bp$V5))
# if any start positions < 0 set to 0
for (i in seq(1,nrow(top_qSNPsgtfdat30bp))){
  if (top_qSNPsgtfdat30bp$V4[i] < 0){
    top_qSNPsgtfdat30bp$V4[i] <- 0
  }
  if (top_qSNPsgtfdat30bp$V5[i] > chrom_dat$length[
    chrom_dat$Chr == as.character(top_qSNPsgtfdat30bp$V1[i])]){
    top_qSNPsgtfdat30bp$V5[i] <- chrom_dat$length[
      chrom_dat$Chr == as.character(top_qSNPsgtfdat30bp$V1[i])]
  }
}
nrow(top_qSNPsgtfdat30bp)
top_qSNPsgtfdat30bp$V4 <- as.integer(as.character(top_qSNPsgtfdat30bp$V4))
top_qSNPsgtfdat30bp$V5 <- as.integer(as.character(top_qSNPsgtfdat30bp$V5))

# ####

#
# Write gtf files of SNPs for other analyses
# ####
write.table(gtfdat, 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_all_snps.gtf",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(gtfdat30bp, 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_all_snps30bp.gtf",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(top_bonfSNPsgtfdat, 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_bonf_snps.gtf",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(top_qSNPsgtfdat, 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_q-value_snps.gtf",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(top_bonfSNPsgtfdat30bp, 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_bonf_snps30bp.gtf",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(top_qSNPsgtfdat30bp, 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_q-value_snps30bp.gtf",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Write SNPs in GOwinda format
head(top_qSNPsgtfdat[,c(1,4)])
nrow(top_qSNPsgtfdat[,c(1,4)])
head(top_bonfSNPsgtfdat[,c(1,4)])
nrow(top_bonfSNPsgtfdat[,c(1,4)])

write.table(top_bonfSNPsgtfdat[,c(1,4)], 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_bonf_snps.tab",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(top_qSNPsgtfdat[,c(1,4)], 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_q-value_snps.tab",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(gtfdat[,c(1,4)], 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_all_snps.tab",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Write list of SNPs IDs
top_bonfSNPs_list <- paste("SNP_",top_bonfSNPs$Chr,
                           "_",top_bonfSNPs$Pos, sep="") 
top_qSNPs_list <- paste("SNP_",top_qSNPs$Chr,"_",
                        top_qSNPs$Pos, sep="") 

allSNPs_list <- paste("SNP_",glmdat$Chr,"_",
                      glmdat$Pos, sep="")

write.table(top_bonfSNPs_list, 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_bonf_snps.list",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(top_qSNPs_list, 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_q-value_snps.list",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(allSNPs_list, 
            file = "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_all_snps.list",
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
# ####





# ####

#
# Check the plots for bonf SNPs
# ####
glmdat_bonf

# ####


