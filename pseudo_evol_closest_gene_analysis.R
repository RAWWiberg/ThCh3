##
# D.pseudoobscura project evolution linse
# Closest Genes
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
options(scipen=6)
# ####
setwd("~/Desktop/Data/RData/pseudo_evol/")
#


# Load Dpse annotation
dpse_ann<-read.table("~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/Dpse_genes_dmelnames.gtf",header = FALSE,sep="\t")
colnames(dpse_ann)<-c("chr","source","type","start","end",
                      "-","-","-","gene")
head(dpse_ann)
# Load the lists of SNPs
top_bonfSNPs_list<-read.table(
  "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_bonf_snps.list",
  sep = "\t", 
  header=FALSE)

top_qSNPs_list<-read.table(
  "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_q-value_snps.list",
  sep = "\t", 
  header=FALSE)

allSNPs_list<-read.table(
  "/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_all_snps.list",
  sep = "\t", 
  header=FALSE)


#
#
#------------------------------------------#
# Look at distribution of SNPs in relation #
# to annotated genes.                      #
#------------------------------------------#
# Load Dpse annotation vs SNP gtf overlaps files
# ####
top_bonf_genes<-read.table("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_bonf_snps_closest_genes_reduced.tab",
                           header = FALSE, sep = "\t")
colnames(top_bonf_genes) <- c("SNP","gene_id","dist")
top_bonf_genes <- top_bonf_genes[!(top_bonf_genes$gene_id == "."),]
nrow(top_bonf_genes)

top_q_genes<-read.table("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_q-value_snps_closest_genes_reduced.tab",
                        header = FALSE, sep = "\t")
colnames(top_q_genes) <- c("SNP","gene_id","dist")
top_q_genes <- top_q_genes[!(top_q_genes$gene_id == "."),]
nrow(top_q_genes)

all_genes <-read.table("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_all_closest_genes_reduced.tab",
                       header = FALSE, sep = "\t")
colnames(all_genes) <- c("SNP","gene_id","dist")
all_genes <- all_genes[
  !(all_genes$gene_id == "." | all_genes$gene_id == "-1"),]
nrow(all_genes)

# distribution of "absolute" distances from gene
top_bonf_dist_hist <- ggplot()+
  geom_histogram(data = top_bonf_genes, aes(abs(dist)/1000))+
  xlab("|Distance from Gene| (kb)")+
  ylab("Count")
top_bonf_dist_hist + my.theme

top_q_dist_hist <- ggplot()+
  geom_histogram(data = top_q_genes, aes(abs(dist)/1000))+
  xlab("|Distance from Gene| (kb)")+
  ylab("Count")
top_q_dist_hist + my.theme

all_dist_hist <- ggplot()+
  geom_histogram(data = all_genes, aes(abs(dist)/1000))+
  xlab("|Distance from Gene| (kb)")+
  ylab("Count")
all_dist_hist + my.theme
# ####


# How many SNPs outside a gene
nrow(top_bonf_genes[abs(top_bonf_genes$dist) > 0,])
nrow(top_bonf_genes[abs(top_bonf_genes$dist) > 0,])/nrow(top_bonf_genes)

nrow(top_q_genes[abs(top_q_genes$dist) > 0,])
nrow(top_q_genes[abs(top_q_genes$dist) > 0,])/nrow(top5_genes)

nrow(all_genes[abs(all_genes$dist) > 0,])
nrow(all_genes[abs(all_genes$dist) > 0,])/nrow(all_genes)

# How many SNPs WITHIN coding region of a gene
nrow(top_bonf_genes[top_bonf_genes$dist == 0,])
nrow(top_bonf_genes[top_bonf_genes$dist == 0,])/nrow(top_bonf_genes)

nrow(top_q_genes[top_q_genes$dist == 0,])
nrow(top_q_genes[top_q_genes$dist == 0,])/nrow(top_q_genes)

nrow(all_genes[all_genes$dist == 0,])
nrow(all_genes[all_genes$dist == 0,])/nrow(all_genes)
# subset to include only genes with SNP within a coding region
top_bonf_genes_withingene <- top_bonf_genes[
  abs(top_bonf_genes$dist) == 0,]

top_q_genes_withingene <- top_q_genes[
  abs(top_q_genes$dist) == 0,]

all_genes_withingene <- all_genes[
  abs(all_genes$dist) == 0,]
# get unique genes
length(unique(top_bonf_genes_withingene$gene_id))

length(unique(top_q_genes_withingene$gene_id))

length(unique(all_genes_withingene$gene_id))


#------------------------------------------------------------#
# How many SNPs within XX (1kb or 1Mb) of a gene
#------------------------------------------------------------#
XX<-1000000
#------------------------------------------------------------#
# UP or DOWNSTREAM, 
# *includes* SNPs that occur WITHIN a gene
#------------------------------------------------------------#
nrow(top_q_genes[
  abs(top_q_genes$dist) <= XX,])
nrow(top_q_genes[
  abs(top_q_genes$dist) <= XX,])/nrow(top_q_genes)

nrow(all_genes[abs(all_genes$dist) <= XX,])
nrow(all_genes[abs(all_genes$dist) <= XX,])/nrow(all_genes)

# Subset to include only genes with SNP within XX
# of a SNP *including* SNPs within gene region
top_q_genes_XX <- top_q_genes[
  abs(top_q_genes$dist) <= XX,]

all_genes_XX <- all_genes[
  abs(all_genes$dist) <= XX,]

length(unique(top_bonf_genes_XX$gene_id))

length(unique(top_q_genes_XX$gene_id))

length(unique(all_genes_XX$gene_id))

# write Gene IDs as lists
write.table(unique(top_bonf_genes_XX$gene_id),
            paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_bonf_closest_genes_",
                  XX,".list",
                  sep = ""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)

write.table(unique(top_q_genes_XX$gene_id),
            paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_q-value_closest_genes_",
                  XX,".list",
                  sep = ""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)

write.table(unique(all_genes_XX$gene_id),
            paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_all_closest_genes_",
                  XX,".list",
                  sep = ""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)


#------------------------------------------------------------#
# UP or DOWNSTREAM, 
# *excludes* SNPs that occur WITHIN a gene
#------------------------------------------------------------#
nrow(top_bonf_genes[
  abs(top_bonf_genes$dist) <= XX & 
    top_bonf_genes$dist != 0,])
nrow(top_bonf_genes[
  abs(top_bonf_genes$dist) <= XX &
    top_bonf_genes$dist != 0,])/nrow(top_bonf_genes)

nrow(top_q_genes[
  abs(top_q_genes$dist) <= XX & 
    top_q_genes$dist != 0,])
nrow(top_q_genes[
  abs(top_q_genes$dist) <= XX &
    top_q_genes$dist != 0,])/nrow(top_q_genes)

nrow(all_genes[
  abs(all_genes$dist) <= XX &
    all_genes$dist != 0,])
nrow(all_genes[
  abs(all_genes$dist) <= XX &
    all_genes$dist != 0,])/nrow(all_genes)

# Subset to include only genes with SNP within XX
# of a SNP *and* not within gene region.
top_bonf_genes_XX_out <- top_bonf_genes[
  abs(top_bonf_genes$dist) <= XX & 
    top_bonf_genes$dist != 0,]

top_q_genes_XX_out <- top_q_genes[
  abs(top_q_genes$dist) <= XX & 
    top_q_genes$dist != 0,]

all_genes_XX_out <- all_genes[
  abs(all_genes$dist) <= XX & 
    all_genes$dist != 0,]

length(unique(top_bonf_genes_XX_out$gene_id))

length(unique(top_q_genes_XX_out$gene_id))

length(unique(all_genes_XX_out$gene_id))

# Write Gene IDs as lists
write.table(unique(top_bonf_genes_XX_out$gene_id),
            paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_bonf_closest_genes_",
                  XX,"_out.list",
                  sep=""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)

write.table(unique(top_q_genes_XX_out$gene_id),
            paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_q-value_closest_genes_",
                  XX,"_out.list",
                  sep=""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)

write.table(unique(all_genes_XX_out$gene_id),
            paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_all_closest_genes_",
                  XX,"_out.list",
                  sep=""),
            quote = FALSE,col.names=FALSE,row.names=FALSE)


#------------------------------------------------------------#
# How many SNPs within XX UPSTREAM of a gene
# *Excludes* SNPs that occur WITHIN a gene
#------------------------------------------------------------#
nrow(top_bonf_genes[top_bonf_genes$dist > -XX & 
                      top_bonf_genes$dist < 0,])
nrow(top_bonf_genes[top_bonf_genes$dist > -XX & 
                      top_bonf_genes$dist < 0,])/nrow(top_bonf_genes)

nrow(top_q_genes[top_q_genes$dist > -XX & 
                   top_q_genes$dist < 0,])
nrow(top_q_genes[top_q_genes$dist > -XX & 
                   top_q_genes$dist < 0,])/nrow(top_q_genes)

nrow(all_genes[all_genes$dist > -XX & 
                 all_genes$dist < 0,])
nrow(all_genes[all_genes$dist > -XX & 
                 all_genes$dist < 0,])/nrow(all_genes)

nrow(top_bonf_genes[top_bonf_genes$dist > -1000,])
nrow(top_bonf_genes[top_bonf_genes$dist > -1000,])/nrow(top_bonf_genes)
# Subset to include only genes with SNP within XX
# of a SNP
top_bonf_genes_XX_ups <- top_bonf_genes[
  top_bonf_genes$dist > -XX & 
    top_bonf_genes$dist < 0,]

top_q_genes_XX_ups <- top_q_genes[
  top_q_genes$dist > -XX & 
    top_q_genes$dist < 0,]


#------------------------------------------------------------#
# Write SNPs within XX UP or DOWNSTREAM of a gene as gtf files
# EXLUDING genes WITHIN the coding region of a gene
#------------------------------------------------------------#
# Load glmdat 
glmdat <- read.table("pseudo_evol_qbglm-unp_cov17-49.rout",
                     sep = "", header = FALSE)
colnames(glmdat) <- c("Chr","Pos","ref",
                      "R1M","R1P","R2M","R2P",
                      "R3M","R3P","R4M","R4P",
                      "t_p")
# t = Treatment (E vs M)[Line]
glmdat$SNP <- paste("SNP_",glmdat$Chr,"_",glmdat$Pos, sep="")
glmdat$Chr3 <- gsub("^", "chrom", glmdat$Chr)
glmdat$Chr2 <- factor(gsub("_group",".",glmdat$Chr))

# Which top SNPs to use
# UP or DOWNSTREAM: EXCLUDING SNPs WITHIN a gene
top_qSNPs<-top_q_genes_XX_out
top_bonfSNPs<-top_bonf_genes_XX_out
dat<-"_out"

# UP or DOWNSTREAM: INCLUDING SNPs WITHIN a gene
top_qSNPs<-top_q_genes_XX
top_bonfSNPs<-top_bonf_genes_XX
dat<-""

# only UPSTREAM: EXCLUDING SNPs WITHIN a gene
top_qSNPs<-top_q_genes_XX_ups
top_bonfSNPs<-top_bonf_genes_XX_ups
dat<-"_ups"

# Subset glmdat to top bonf SNPs 
glmdat_bonf<-glmdat[glmdat$SNP %in% top_bonfSNPs$SNP,]
head(glmdat_bonf)

length(top_bonfSNPs$SNP)
nrow(glmdat_bonf)
# Make a gtf file
top_bonf_SNPsgtfdat30bp <- cbind(as.character(glmdat_bonf$Chr),
                                 rep("PoPoolation2",nrow(glmdat_bonf)),
                                 rep("CDS",nrow(glmdat_bonf)),
                                 (glmdat_bonf$Pos-30),
                                 (glmdat_bonf$Pos+30),
                                 rep(".",nrow(glmdat_bonf)),
                                 rep(".",nrow(glmdat_bonf)),
                                 rep(".",nrow(glmdat_bonf)),
                                 paste(
                                   paste('gene_id "SNP_',
                                         seq(1,nrow(glmdat_bonf),1),
                                         '"',sep=""),
                                   paste('transcript_id "SNP_',
                                         seq(1,nrow(glmdat_bonf),1),
                                         '"',sep=""),
                                   paste('T-p_value "',
                                         glmdat_bonf$t_p,'"',sep=""),
                                   sep = "; "))
top_bonf_SNPsgtfdat30bp<-as.data.frame(top_bonf_SNPsgtfdat30bp)
top_bonf_SNPsgtfdat30bp$V4 <- as.integer(
  as.character(top_bonf_SNPsgtfdat30bp$V4)) 
top_bonf_SNPsgtfdat30bp$V5 <- as.integer(
  as.character(top_bonf_SNPsgtfdat30bp$V5)) 

write.table(top_bonf_SNPsgtfdat30bp, 
            file = paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_bonf_snps30bp_",XX,
                         dat,".gtf",
                         sep=""),
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Subset glmdat to top bonf SNPs
glmdat_q<-glmdat[glmdat$SNP %in% top_qSNPs$SNP,]
head(glmdat_q)
length(top_qSNPs)
nrow(glmdat_q)
# Make a gtf file
top_q_SNPsgtfdat30bp <- cbind(as.character(glmdat_q$Chr),
                                 rep("PoPoolation2",nrow(glmdat_q)),
                                 rep("CDS",nrow(glmdat_q)),
                                 (glmdat_q$Pos-30),
                                 (glmdat_q$Pos+30),
                                 rep(".",nrow(glmdat_q)),
                                 rep(".",nrow(glmdat_q)),
                                 rep(".",nrow(glmdat_q)),
                                 paste(
                                   paste('gene_id "SNP_',
                                         seq(1,nrow(glmdat_q),1),
                                         '"',sep=""),
                                   paste('transcript_id "SNP_',
                                         seq(1,nrow(glmdat_q),1),
                                         '"',sep=""),
                                   paste('T-p_value "',
                                         glmdat_q$t_p,'"',sep=""),
                                   sep = "; "))
top_q_SNPsgtfdat30bp<-as.data.frame(top_q_SNPsgtfdat30bp)
top_q_SNPsgtfdat30bp$V4 <- as.integer(
  as.character(top_q_SNPsgtfdat30bp$V4)) 
top_q_SNPsgtfdat30bp$V5 <- as.integer(
  as.character(top_q_SNPsgtfdat30bp$V5)) 

write.table(top_q_SNPsgtfdat30bp, 
            file = paste("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/pseudo_evol_glm_top_q_snps30bp_",XX,
                         dat,".gtf",
                         sep=""),
            sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


#------------------------------------#
# Are the genes near top SNPs among:
#   1) Candidate genes
#   2) DE genes
#   3) Peak genes
#------------------------------------#
top_qSNPs
length(unique(top_qSNPs$gene_id))
top_bonfSNPs
length(unique(top_bonfSNPs$gene_id))
XX
dat

#----------------#
# CANDIDATE GENES
#----------------#
# Load the candidate genes
cand_genes <- read.table(
  "~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/TableS2_pseudoobscura_remating_candidate_genes.csv",
  header = TRUE,sep=",")
head(cand_genes)
cand_genes_V<-gsub("gn","GN",cand_genes$D..melanogaster.FlyBase.ID)
# Pickpocket genes from Paris
pkk <- c("FBGN0034730","FBGN0034489")

top_bonfSNPs$gene_id<-gsub("_.","",top_bonfSNPs$gene_id)
top_qSNPs$gene_id<-gsub("_.","",top_qSNPs$gene_id)

cand_genes_V[which(cand_genes_V %in% unique(top_bonfSNPs$gene_id))]
pkk[which(pkk %in% unique(top_bonfSNPs$gene_id))]

cand_genes_V[which(cand_genes_V %in% unique(top_qSNPs$gene_id))]
pkk[which(pkk %in% unique(top_qSNPs$gene_id))]

#----------------#
# PEAK GENES
#----------------#
# Load the "peak genes"
peak_genes<-read.table("~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/TableS3_dpse_glm_bonf_peaks_regions.tab",
                       header=TRUE,sep="\t")

# Which peak genes have a SNP nearby
peak_genes_w_q_snp <-peak_genes[
  peak_genes$FBid %in% gsub("FBGN","FBgn",top_qSNPs$gene_id),] 

length(unique(peak_genes_w_bonf_snp$FBid))
length(unique(peak_genes_w_q_snp$FBid))

# Write the genes with nearby snps as a table
write.table(peak_genes_w_q_snp,
            paste("~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/Table3.2_dpse_glm_q_peaks_regions_wSNPs_",
                  XX,
                  dat,".tab",sep=""),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

#----------------#
# NOT PEAK GENES
#----------------#
# Which genes are OUTSIDE the peak regions but within XX b of
# a gene.
not_peak_genes_w_q_snp<-top_qSNPs$gene_id[
  which(
    !(gsub("FBGN","FBgn",top_qSNPs$gene_id) %in% 
        peak_genes_w_q_snp$FBid))]

length(unique(not_peak_genes_w_q_snp))
dpse_ann$gene<-gsub("gene_id ","",dpse_ann$gene)
dpse_ann$gene<-gsub(";","",dpse_ann$gene)
dpse_ann$gene<-gsub("_.","",dpse_ann$gene)

not_peak_genes_w_q_snp_tab<-dpse_ann[
  dpse_ann$gene %in% not_peak_genes_w_q_snp,
  c(1,4,5,9)]
colnames(not_peak_genes_w_q_snp_tab)<-c("chr","start","end","gene")
not_peak_genes_w_q_snp_tab$gene<-gsub("FBGN","FBgn",
                                      not_peak_genes_w_q_snp_tab$gene)

# Write genes with nearby snps as list.
write.table(not_peak_genes_w_q_snp_tab$gene,
            paste("~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/dpse_glm_q_not_peak_genes_wSNPs",
                  XX,dat,".list",sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

# Submit the above lists to flybase to get the actual gene name
# then load the conversion table to R and add the actual gene name
# to tables S4.1 and S4.2
not_peak_genes_w_q_snp_tab_conv<-read.table("~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/dpse_glm_q_not_peak_genes_wSNPs1000000_FBconv.tab",
                          header=FALSE,sep="\t")
colnames(not_peak_genes_w_q_snp_tab_conv)<-c("geneidsubmit","FBgeneid",
                           "FBgeneidconv","genename")
not_peak_genes_w_q_snp_tab_conv

# Add the converted gene name to the table
not_peak_genes_w_q_snp_tab$name<-vector(
  length=nrow(not_peak_genes_w_q_snp_tab))
for(id in not_peak_genes_w_q_snp_tab$gene){
  genename<-as.character(not_peak_genes_w_q_snp_tab_conv$genename[
    not_peak_genes_w_q_snp_tab_conv$geneidsubmit == id])
  not_peak_genes_w_q_snp_tab$name[
    not_peak_genes_w_q_snp_tab$gene == id]<-genename  
}
not_peak_genes_w_q_snp_tab

# Write genes as a table
write.table(not_peak_genes_w_q_snp_tab,
            paste("~/Desktop/Data/pseudo_project/evol_lines/genomics/SNP_analysis/TableS4.2_dpse_glm_q_not_peak_genes_wSNPs",
                  XX,dat,".tab",sep=""),
            quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

#----------------#
# DE GENES
#----------------#
peak_genes_w_q_snp<-peak_genes_w_q_snp[c(1,2,3,4,5)]
peak_genes_w_q_snp
not_peak_genes_w_q_snp_tab

all_genes_w_qsnps<-c(as.character(peak_genes_w_q_snp$name),
                    as.character(not_peak_genes_w_q_snp_tab$name))
length(all_genes_w_qsnps)

dmel_table <- read.table("/home/axel/Desktop/Data/pseudo_project/reference/dpse_r31_FB2013_02/gff/dpse_FlyBase_orth_dmel_id.tab",
                         sep = "\t",
                         quote = "\"",
                         header = FALSE)
colnames(dmel_table) <- c("name","FB_id1", "FB_id2","name2")
head(dmel_table)
dpse_table <- read.table("/home/axel/Desktop/Data/pseudo_project/reference/dpse_r31_FB2013_02/gff/dpse_FlyBase_Fields_download.txt",
                         sep = "\t",
                         quote = "\"",
                         header = FALSE)
colnames(dpse_table) <- c("FB_id1", "FB_id2","name","name2","species")
head(dpse_table)
nrow(dpse_table)

all_genes_w_qsnps_dpse_names<-dpse_table$FB_id1[
  which(dpse_table$name %in% all_genes_w_qsnps)]

#
# IMMONEN ET AL., 2014
#
# Read in Immonen et al., 2014 data
immonen<-read.table("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/EI_gene_expression/Immonen_et_al_2014_Supp_Tables_diff_exprs_genes_EM.csv",
                    header=TRUE,sep=",")
colnames(immonen)<-c(colnames(immonen)[1:6],"FBid")
head(immonen)

# How many of Immonen et al., 2014 DE genes are 
# among the all_genes_w_snps
# let's call this "enrichment score"(ES)
length(immonen$FBid)
immonen_ES<-length(immonen$FBid[
  which(immonen$FBid %in% all_genes_w_qsnps_dpse_names)])
immonen_ES
immonen_ES/length(all_genes_w_qsnps_dpse_names)
# Which genes
immonen_genes<-immonen$FBid[
  which(immonen$FBid %in% all_genes_w_qsnps_dpse_names)]
dpse_table[which(dpse_table$FB_id1 %in% immonen_genes),]
# How many aren't
length(immonen$FBid)-length(immonen$FBid[
  which(immonen$FBid %in% all_genes_w_qsnps_dpse_names)])

#
# VELTSOS ET AL., UNPUBLISHE
#
# Read in Veltsos data
veltsos<-read.table("/home/axel/Desktop/Data/pseudo_project/evol_lines/genomics/PV_gene_expression/PV_Contrasts_MFBHVC_all_DE_genes.tab",
                    header=FALSE,sep=" ")
colnames(veltsos)<-c("FBid","logFC","p-value","FDR","name","contrast")
head(veltsos)
str(veltsos)
nrow(veltsos)
hist(veltsos$FDR)
hist(veltsos$logFC)
# Subset to a particular FDR
FDR<-0.1
veltsos_FDR<-veltsos[veltsos$FDR<FDR,]

# How many of Veltsos DE genes are among the all_genes_w_snps
# let's call this "enrichment score"(ES)
length(veltsos$FBid)
length(veltsos_FDR$FBid)
length(all_genes_w_qsnps_dpse_names)
veltsos_ES<-length(veltsos_FDR$FBid[
  which(veltsos_FDR$FBid %in% all_genes_w_qsnps_dpse_names)])
veltsos_ES
veltsos_ES/length(all_genes_w_qsnps_dpse_names)
# Which genes?
veltsos_genes<-veltsos_FDR[
  which(veltsos_FDR$FBid %in% all_genes_w_qsnps_dpse_names),]
# How many aren't
length(veltsos_FDR$FBid)-length(veltsos_FDR$FBid[
  which(veltsos_FDR$FBid %in% all_genes_w_qsnps_dpse_names)])





# RUN A BOOTSTRAPPING ROUTINE
boots<-10000
# Sample an equal number of genes randomly from the dpse_table
N<-length(all_genes_w_qsnps_dpse_names)
dataset<-veltsos_FDR
ES<-veltsos_ES
# How many of Random genes are among the all_genes_w_snps
boot_samps<-vector(length=boots)
b<-0
while(b <= boots){
  randm<-dpse_table$FB_id1[
    sample(seq(1,nrow(dpse_table),1),N,replace = FALSE)]
  boot_samps[b]<-length(dataset$FBid[
    which(dataset$FBid %in% randm)])
  b<-b+1
}

# What proportion of the boot_samps are >= ES
length(boot_samps[boot_samps>=ES])/length(boot_samps)
mean(boot_samps)/length(all_genes_w_qsnps_dpse_names)
hist(boot_samps)

# CHECK OVERLAP BETWEEN THE DATASETS

# How many in common between immonen and veltsos
# (veltsos$FBid is the shorter vector)
length(veltsos$FBid[which(veltsos_FDR$FBid %in% immonen$FBid)])
length(unique(veltsos$FBid[which(veltsos_FDR$FBid %in% immonen$FBid)]))
127/length(unique(veltsos_FDR$FBid))
127/length(unique(immonen$FBid))

# Are the immonen_genes the same as the veltsos_genes?
# (veltsos_genes is the shorter vector)
veltsos_genes[which(veltsos_genes$FBid %in% immonen_genes),]






#
# ADDITIONAL ANALYSIS IN VELTSOS ET AL., UNPUBLISHED
#
# Get the p-values and q-values for each of the SNPs
top_q_genes_1000000bp_pvals<-glmdat_q[
  which(glmdat_q$SNP %in% top_q_genes_1000000bp$SNP),
  c(13,12,16)]

# Make sure the SNPs are in the same order
top_q_genes_1000000bp_pvals<-top_q_genes_1000000bp_pvals[
  order(top_q_genes_1000000bp_pvals$SNP),]
top_q_genes_1000000bp<-top_q_genes_1000000bp[
  order(top_q_genes_1000000bp$SNP),]
# Get the mean p-value for each gene
top_q_genes_1000000bp<-cbind(top_q_genes_1000000bp,
                             top_q_genes_1000000bp_pvals[,c(2,3)])
str(top_q_genes_1000000bp)
mean_p_gene<-aggregate(top_q_genes_1000000bp$t_p,
                       by = list(top_q_genes_1000000bp$gene_id),
                       FUN=mean)
colnames(mean_p_gene)<-c("gene_id","t_p")
mean_p_gene$gene_id<-gsub("_.","",mean_p_gene$gene_id)
mean_p_gene$gene_id<-gsub("FBGN","FBgn",mean_p_gene$gene_id)
# For lack of a better solution I exclude duplicates
duplicates<-mean_p_gene$gene_id[duplicated(mean_p_gene$gene_id)]
mean_p_gene<-mean_p_gene[which(!(mean_p_gene$gene_id %in% duplicates)),]
head(mean_p_gene)
mean_p_gene

# Get the max p-value for each gene
min_p_gene<-aggregate(top_q_genes_1000000bp$t_p,
                      by = list(top_q_genes_1000000bp$gene_id),
                      FUN=min)
colnames(min_p_gene)<-c("gene_id","t_p")
min_p_gene$gene_id<-gsub("_.","",min_p_gene$gene_id)
min_p_gene$gene_id<-gsub("FBGN","FBgn",min_p_gene$gene_id)
# For lack of a better solution I exclude duplicates
duplicates<-min_p_gene$gene_id[duplicated(min_p_gene$gene_id)]
min_p_gene<-min_p_gene[which(!(min_p_gene$gene_id %in% duplicates)),]

# Check the correlation between the mean and the max
cor.test(min_p_gene$t_p,mean_p_gene$t_p,method = "spearman")
cor.test(min_p_gene$t_p,mean_p_gene$t_p,method = "pearson")
plot(min_p_gene$t_p,mean_p_gene$t_p)

# Plot QBGLM p-values against DE p-values
# Get gene names
mean_p_gene$name<-vector(length=nrow(mean_p_gene))
mean_p_gene$dpse_gene_id<-vector(length=nrow(mean_p_gene))
for(gene in mean_p_gene$gene_id){
  name<-as.character(
    dmel_table$name[dmel_table$FB_id1 == gsub("FBgn","FBGN",gene)])
  mean_p_gene$name[mean_p_gene$gene_id == gene]<-name
  dpse_name<-as.character(dpse_table$FB_id1[
    dpse_table$name == name])[1]
  mean_p_gene$dpse_gene_id[mean_p_gene$gene_id == gene]<-dpse_name
}
mean_p_gene
# Get the DE p-values for SNP associated genes
# Make collection data.frame
out_data<-data.frame()
for(contrast in unique(veltsos$contrast)){
  print(contrast)
  veltsos_sub<-veltsos[veltsos$contrast==contrast,]
  # Subset veltsos sub to get subset data.
  veltsos_sub<-veltsos_sub[
    which(veltsos_sub$FBid %in% mean_p_gene$dpse_gene_id),]
  # Order by gene name
  veltsos_sub<-veltsos_sub[order(veltsos_sub$FBid),]
  # Subset mean_p_gene to only those in the veltsos data
  out_mean_p_gene<-mean_p_gene[
    which(mean_p_gene$dpse_gene_id %in% veltsos_sub$FBid),]
  # Order by gene name
  out_mean_p_gene<-out_mean_p_gene[order(out_mean_p_gene$dpse_gene_id),]
  # Collect data
  coldat<-cbind(out_mean_p_gene,veltsos_sub[,3],
                rep(contrast,nrow(out_mean_p_gene)))
  # Save data
  out_data<-rbind(out_data,coldat)
}
colnames(out_data)<-c(colnames(mean_p_gene)[1:4],"DE_pval","contrast")

head(out_data)
ggplot()+
  geom_point(data=out_data,aes(x=-log10(t_p),y=-log10(DE_pval),
                               group=contrast))+
  geom_hline(yintercept = -log10(0.05),linetype="dashed",colour="red")+
  geom_vline(xintercept = -log10(0.05),linetype="dashed",colour="red")+
  xlim(0,10)+
  ylim(0,10)+
  facet_wrap(~contrast)
tapply(out_data$gene_id,INDEX = list(out_data$contrast),length)
length(unique(out_data$gene_id))

write.table(out_data[out_data$DE_pval < 0.05,],
            "~/Desktop/PV_AW_comparison.txt",row.names = FALSE,
            quote=FALSE)
