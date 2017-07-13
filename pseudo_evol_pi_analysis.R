##
# D.pseudoobscura project evolution linse
# PoPoolation pi data analysis
# Last Modified: March 2017
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
se <- function(x) {sqrt(var(x, na.rm = TRUE))/sqrt(length(x))}
setwd("~/Desktop/Data/RData/pseudo_evol/")
source("~/Desktop/Data/RData/RScripts/ggplot_theme.R")
options(scipen = 6)


# Read in data
# ####
###
# E1
###
E1pi<-read.table("pseudo_evol_E1.pi",header = FALSE,
                  sep = "\t")
colnames(E1pi)<-c("Chr","WinPos","NrSNPs","PropGdReads","pi")
E1pi$pi <- as.character(E1pi$pi)
# Subset to only known chromosomes
E1pi<-E1pi[grep("Unknown",E1pi$Chr,invert = TRUE),]
E1pi$tr<-rep("E",nrow(E1pi))
E1pi$rep<-rep("1",nrow(E1pi))
str(E1pi)
head(E1pi)
nrow(E1pi)
###
# M1
###
M1pi<-read.table("pseudo_evol_M1.pi",header = FALSE,
                 sep = "\t")
colnames(M1pi)<-c("Chr","WinPos","NrSNPs","PropGdReads","pi")
M1pi$pi <- as.character(M1pi$pi)
# Subset to only known chromosomes
M1pi<-M1pi[grep("Unknown",M1pi$Chr,invert = TRUE),]
M1pi$tr<-rep("M",nrow(M1pi))
M1pi$rep<-rep("1",nrow(M1pi))
head(M1pi)
nrow(M1pi)
###
# E2
###
E2pi<-read.table("pseudo_evol_E2.pi",header = FALSE,
                 sep = "\t")
colnames(E2pi)<-c("Chr","WinPos","NrSNPs","PropGdReads","pi")
E2pi$pi <- as.character(E2pi$pi)
# Subset to only known chromosomes
E2pi<-E2pi[grep("Unknown",E2pi$Chr,invert = TRUE),]
E2pi$tr<-rep("E",nrow(E2pi))
E2pi$rep<-rep("2",nrow(E2pi))
str(E2pi)
head(E2pi)
nrow(E2pi)
###
# M2
###
M2pi<-read.table("pseudo_evol_M2.pi",header = FALSE,
                 sep = "\t")
colnames(M2pi)<-c("Chr","WinPos","NrSNPs","PropGdReads","pi")
M2pi$pi <- as.character(M2pi$pi)
# Subset to only known chromosomes
M2pi<-M2pi[grep("Unknown",M2pi$Chr,invert = TRUE),]
M2pi$tr<-rep("M",nrow(M2pi))
M2pi$rep<-rep("2",nrow(M2pi))
head(M2pi)
nrow(M2pi)
###
# E3
###
E3pi<-read.table("pseudo_evol_E3.pi",header = FALSE,
                 sep = "\t")
colnames(E3pi)<-c("Chr","WinPos","NrSNPs","PropGdReads","pi")
E3pi$pi <- as.character(E3pi$pi)
# Subset to only known chromosomes
E3pi<-E3pi[grep("Unknown",E3pi$Chr,invert = TRUE),]
E3pi$tr<-rep("E",nrow(E3pi))
E3pi$rep<-rep("3",nrow(E3pi))
str(E3pi)
head(E3pi)
nrow(E3pi)
###
# M3
###
M3pi<-read.table("pseudo_evol_M3.pi",header = FALSE,
                 sep = "\t")
colnames(M3pi)<-c("Chr","WinPos","NrSNPs","PropGdReads","pi")
M3pi$pi <- as.character(M3pi$pi)
# Subset to only known chromosomes
M3pi<-M3pi[grep("Unknown",M3pi$Chr,invert = TRUE),]
M3pi$tr<-rep("M",nrow(M3pi))
M3pi$rep<-rep("3",nrow(M3pi))
head(M3pi)
nrow(M3pi)
###
# E4
###
E4pi<-read.table("pseudo_evol_E4.pi",header = FALSE,
                 sep = "\t")
colnames(E4pi)<-c("Chr","WinPos","NrSNPs","PropGdReads","pi")
E4pi$pi <- as.character(E4pi$pi)
# Subset to only known chromosomes
E4pi<-E4pi[grep("Unknown",E4pi$Chr,invert = TRUE),]
E4pi$tr<-rep("E",nrow(E4pi))
E4pi$rep<-rep("4",nrow(E4pi))
str(E4pi)
head(E4pi)
nrow(E4pi)
###
# M4
###
M4pi<-read.table("pseudo_evol_M4.pi",header = FALSE,
                 sep = "\t")
colnames(M4pi)<-c("Chr","WinPos","NrSNPs","PropGdReads","pi")
M4pi$pi <- as.character(M4pi$pi)
# Subset to only known chromosomes
M4pi<-M4pi[grep("Unknown",M4pi$Chr,invert = TRUE),]
M4pi$tr<-rep("M",nrow(M4pi))
M4pi$rep<-rep("4",nrow(M4pi))
head(M4pi)
nrow(M4pi)
# ####

# Combine pi sets
# pi for E and M as separate columns
# ####
E_vs_M_pi<-cbind(as.character(c(M1pi$pi,M2pi$pi,M3pi$pi,M4pi$pi)),
                 as.character(c(E1pi$pi,E2pi$pi,E3pi$pi,E4pi$pi)),
                 c(rep("1",length(M1pi$pi)),
                   rep("2",length(M2pi$pi)),
                   rep("3",length(M3pi$pi)),
                   rep("4",length(M4pi$pi))),
                 c(as.character(E1pi$Chr),
                   as.character(E2pi$Chr),
                   as.character(E3pi$Chr),
                   as.character(E4pi$Chr)))

E_vs_M_pi<-as.data.frame(E_vs_M_pi)
colnames(E_vs_M_pi)<-c("Mpi","Epi","Rep","Chr")
E_vs_M_pi<-E_vs_M_pi[E_vs_M_pi$Epi != "na" & E_vs_M_pi$Mpi != "na",]
E_vs_M_pi$Mpi<-as.numeric(as.character(E_vs_M_pi$Mpi))
E_vs_M_pi$Epi<-as.numeric(as.character(E_vs_M_pi$Epi))
E_vs_M_pi$Chr<-gsub("_group",".",E_vs_M_pi$Chr)
head(E_vs_M_pi)
str(E_vs_M_pi)
nrow(E_vs_M_pi)
# ####

# pi for E an M (melted)
# ####
pidat<-rbind(E1pi,M1pi,
             E2pi,M2pi,
             E3pi,M3pi,
             E4pi,M4pi)
# Remove "na"
#pidat<-pidat[pidat$pi != "na",]
pidat$pi<-as.numeric(pidat$pi)
pidat$chrty<-vector(length=nrow(pidat))
for (row in seq(1,nrow(pidat))){
  if (grepl("X",as.character(pidat$Chr[row]))){
    pidat$chrty[row]<-"X"
    }
  else{
    pidat$chrty[row]<-"A"
  }
}

# What is average diversity across all chromosomes
tapply(pidat$pi,INDEX = list(pidat$tr,pidat$chrty,pidat$rep),
       mean,na.rm=TRUE)
# What is the error
tapply(pidat$pi,INDEX = list(pidat$tr,pidat$chrty,pidat$rep),se)*1.69

# E lines seem to have lower diversity overall:
# indicative of selective sweeps (or BG selection)?

Epidat<-rbind(E1pi,
              E2pi,
              E3pi,
              E4pi)
Epidat<-Epidat[Epidat$pi != "na",]
Epidat$pi<-as.numeric(Epidat$pi)
Epidat$Chr<-factor(Epidat$Chr)
Mpidat<-rbind(M1pi,
              M2pi,
              M3pi,
              M4pi)
Mpidat<-Mpidat[Mpidat$pi != "na",]
Mpidat$pi<-as.numeric(Mpidat$pi)
Mpidat$Chr<-factor(Mpidat$Chr)
head(Mpidat)
# Get mean pi across treatments
###
# E
###
Emeandata<-tapply(Epidat$pi,
                 INDEX = list(Epidat$WinPos,Epidat$Chr),mean,na.rm=TRUE)
head(Emeandata)
winpos<-dimnames(Emeandata)[[1]]
chrs<-dimnames(Emeandata)[[2]]
Emeandata<-as.data.frame(Emeandata)
Emeandata$WinPos<-rownames(Emeandata)
head(Emeandata)
str(Emeandata)

Ecombmeandata<-data.frame(winpos=vector(),pi=vector,Chr=vector())
for(chr in chrs){
  pi<-Emeandata[,chr]
  chr<-rep(chr,length(pi))
  winpos<-as.numeric(rownames(Emeandata))
  data<-data.frame(winpos=winpos,pi=pi,Chr=chr)
  Ecombmeandata<-rbind(Ecombmeandata,data)
}
Ecombmeandata<-na.omit(Ecombmeandata)
Ecombmeandata$pi<-as.numeric(Ecombmeandata$pi)
rownames(Ecombmeandata)<-seq(1,nrow(Ecombmeandata))
Ecombmeandata$tr<-rep("E",nrow(Ecombmeandata))
str(Ecombmeandata)
head(Ecombmeandata)
unique(Ecombmeandata$Chr)

###
# M
###
Mmeandata<-tapply(Mpidat$pi,
                 INDEX = list(Mpidat$WinPos,Mpidat$Chr),mean,na.rm=TRUE)
str(Mmeandata)
winpos<-dimnames(Mmeandata)[[1]]
chrs<-dimnames(Mmeandata)[[2]]
Mmeandata<-as.data.frame(Mmeandata)
Mmeandata$WinPos<-rownames(Mmeandata)
head(Mmeandata)
str(Mmeandata)

Mcombmeandata<-data.frame(winpos=vector(),pi=vector,Chr=vector())
for(chr in chrs){
  pi<-Mmeandata[,chr]
  chr<-rep(chr,length(pi))
  winpos<-as.numeric(rownames(Mmeandata))
  data<-data.frame(winpos=winpos,pi=pi,Chr=chr)
  Mcombmeandata<-rbind(Mcombmeandata,data)
}
Mcombmeandata<-na.omit(Mcombmeandata)
Mcombmeandata$pi<-as.numeric(Mcombmeandata$pi)
rownames(Mcombmeandata)<-seq(1,nrow(Mcombmeandata))
Mcombmeandata$tr<-rep("M",nrow(Mcombmeandata))
head(Mcombmeandata)

unique(Mcombmeandata$Chr)

pimeandat<-rbind(Mcombmeandata,Ecombmeandata)
pimeandat$Chr<-gsub("_group",".",pimeandat$Chr)
head(pimeandat)
head(pidat)
pidat$Chr<-gsub("_group",".",pidat$Chr)
# ####

# Plot
# Colours
pal <- c('#d7191c','#fdae61','#abd9e9','#2c7bb6')

# Add chimney regions
chim_reg_dat<-data.frame(Chr=c("2",
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
# Read in the top SNPs
glmdat_top_q <- read.table(
  "pseudo_evol_qbglm-unp_minc16_cov17-49_top_qval_snps.rout",
                     sep = "", header = TRUE)
head(glmdat_top_q)

# Plot
ggplot()+
  geom_line(data=pidat,
            aes(WinPos/1000000,pi,colour=rep),alpha=1/2)+
  geom_line(data=pimeandat,aes(winpos/1000000,pi),size=0.5)+
  xlab("")+
  ylab(expression(pi))+
  scale_y_continuous(limits=c(0,0.035),breaks = c(0,0.015,0.03))+
  scale_colour_manual("Replicate",
                      values=pal)+
  scale_x_continuous(breaks=c(0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30),
                     labels=c("0","","5","","10","","15","","20","","25","","30"))+
  facet_grid(Chr~tr)+
  my.theme+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=11),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "top",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white",colour="grey"))


# Ratios of pi on X versus Autosomes
head(pidat)
A_pidat<-pidat[grep("X",pidat$Chr,invert=TRUE),]
A_pidat$Chr<-gsub("\\..*","",A_pidat$Chr)
A_chr_m<-tapply(A_pidat$pi,
                INDEX=list(A_pidat$tr,A_pidat$rep,A_pidat$Chr),mean,na.rm=TRUE)
A_chr_m_d<-data.frame(pi=matrix(A_chr_m,ncol = 1,byrow = TRUE),
                      tr=dimnames(A_chr_m)[[1]],
                      rep=rep(as.character(dimnames(A_chr_m)[[2]]),each=2),
                      chr=rep(as.character(dimnames(A_chr_m)[[3]]),each=4))
X_pidat<-pidat[grep("X",pidat$Chr),]
X_pidat$Chr<-gsub(".\\..*","",X_pidat$Chr)
X_chr_m<-tapply(X_pidat$pi,
                INDEX=list(X_pidat$tr,X_pidat$rep,X_pidat$Chr),mean,na.rm=TRUE)
X_chr_m_d<-data.frame(pi=matrix(X_chr_m,ncol = 1,byrow = TRUE),
                      tr=dimnames(X_chr_m)[[1]],
                      rep=rep(as.character(dimnames(X_chr_m)[[2]]),each=2),
                      chr=rep(as.character(dimnames(X_chr_m)[[3]]),each=4))
Api<-tapply(A_chr_m_d$pi,INDEX=list(A_chr_m_d$tr,A_chr_m_d$rep),mean)
Apisd<-tapply(A_chr_m_d$pi,INDEX=list(A_chr_m_d$tr,A_chr_m_d$rep),sd)

Xpi<-tapply(X_chr_m_d$pi,INDEX=list(X_chr_m_d$tr,X_chr_m_d$rep),mean)
Xpisd<-tapply(X_chr_m_d$pi,INDEX=list(X_chr_m_d$tr,X_chr_m_d$rep),sd)

Api
Apid<-data.frame(pi=matrix(Api,ncol=1),
                 pisd=matrix(Apisd,ncol=1),
                 rep=rep(dimnames(Api)[[2]],each=2),
           tr=rep(dimnames(Api)[[1]]))
Apid
Xpi
Xpid<-data.frame(pi=matrix(Xpi,ncol=1),
                 pisd=matrix(Xpisd,ncol=1),
                 rep=rep(dimnames(Xpi)[[2]],each=2),
                 tr=rep(dimnames(Xpi)[[1]]))
Xpid
# AUTOSOMES
ggplot()+
  geom_bar(data=Apid,aes(rep,pi,fill=tr),
           position="dodge",
           stat="identity")+
  geom_errorbar(data=Apid,
                aes(rep,ymin=pi-pisd,ymax=pi+pisd,group=tr),
                width=0,position=position_dodge(0.9))+
  xlab("Replicate")+
  ylab(expression(pi))+
  my.theme+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=11),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "right",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white",colour="grey"))

t.test(Apid$pi~Apid$tr)
# Exclude rep 2
t.test(Apid$pi[Apid$rep!="2"]~Apid$tr[Apid$rep!="2"])
ovApi<-tapply(Apid$pi,INDEX=list(Apid$tr),mean)
ovApise<-2*tapply(Apid$pi,INDEX=list(Apid$tr),se)
ovApid<-data.frame(pi=matrix(ovApi,ncol=1),se=matrix(ovApise,ncol=1),
                  tr=dimnames(ovApi)[[1]])
ggplot()+
  geom_point(data=ovApid,aes(tr,pi),
           size=3)+
  geom_errorbar(data=ovApid,
                aes(tr,ymin=pi-se,ymax=pi+se,group=tr),
                width=0,position=position_dodge(0.9))+
  my.theme+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=11),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "right",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white",colour="grey"))

# XSOMES
ggplot()+
  geom_bar(data=Xpid,aes(rep,pi,fill=tr),
           position="dodge",
           stat="identity")+
  xlab("Replicate")+
  ylab(expression(pi))+
  my.theme+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=11),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "right",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white",colour="grey"))

t.test(Xpid$pi~Xpid$tr)
# Exclude rep 2
t.test(Xpid$pi[Xpid$rep!="2"]~Xpid$tr[Xpid$rep!="2"])
ovXpi<-tapply(Xpid$pi,INDEX=list(Xpid$tr),mean)
ovXpise<-2*tapply(Xpid$pi,INDEX=list(Xpid$tr),se)
ovXpid<-data.frame(pi=matrix(ovXpi,ncol=1),se=matrix(ovXpise,ncol=1),
                   tr=dimnames(ovXpi)[[1]])
ggplot()+
  geom_point(data=ovXpid,aes(tr,pi),
             size=3)+
  geom_errorbar(data=ovXpid,
                aes(tr,ymin=pi-se,ymax=pi+se,group=tr),
                width=0,position=position_dodge(0.9))+
  my.theme+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=11),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "right",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white",colour="grey"))


# RATIO OF X:A pi
ovXpid$pi/ovApid$pi
Apid$Xpi<-Xpid$pi
ggplot()+
  geom_bar(data=Apid,aes(rep,Xpi/pi,fill=tr),
           position="dodge",
           stat="identity")+
  xlab("Replicate")+
  ylab(expression(paste(pi[X],":",pi[A]," Ratio")))+
  my.theme+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size=11),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "right",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white",colour="grey"))

tapply(Apid$Xpi/Apid$pi,INDEX=list(Apid$tr),mean)
2*tapply(Apid$Xpi/Apid$pi,INDEX=list(Apid$tr),se)

tapply(Apid$Xpi[Apid$rep!="2"]/Apid$pi[Apid$rep!="2"],
       INDEX=list(Apid$tr[Apid$rep!="2"]),mean)
2*tapply(Apid$Xpi[Apid$rep!="2"]/Apid$pi[Apid$rep!="2"],
         INDEX=list(Apid$tr[Apid$rep!="2"]),se)

#Extract the data for each chimney region
chim_reg_dat
head(pidat)
head(pimeandat)
head(glmdat_top_q)
glmdat_top_q$Chr2<-gsub("chrom","",glmdat_top_q$Chr2)
zoompidat<-data.frame()
zoompimeandat<-data.frame()
zoomsnpdata<-data.frame()
for(reg in chim_reg_dat$reg){
  chr<-as.character(chim_reg_dat$Chr[chim_reg_dat$reg==reg])
  start<-chim_reg_dat$start[chim_reg_dat$reg==reg]-50000
  end<-chim_reg_dat$end[chim_reg_dat$reg==reg]+50000
  region<-paste(chr,":",start/1000000,"-",end/1000000,sep="")
  
  regiondat<-pidat[pidat$Chr==chr & pidat$WinPos>start & 
                     pidat$WinPos<end,]
  regiondat<-cbind(regiondat,
                   rep(region,nrow(regiondat)))
  
  regionmeandat<-pimeandat[pimeandat$Chr==chr & pimeandat$winpos>start & 
                             pimeandat$winpos<end,]
  regionmeandat<-cbind(regionmeandat,
                       rep(region,nrow(regionmeandat)))
  regionsnpsdat<-glmdat_top_q[glmdat_top_q$Chr2==chr & 
                                glmdat_top_q$Pos > start & 
                                glmdat_top_q$Pos < end,]
  regionsnpsdat<-cbind(regionsnpsdat,
                       rep(region,nrow(regionsnpsdat)))
  zoompidat<-rbind(zoompidat,
                   regiondat)
  zoompimeandat<-rbind(zoompimeandat,
                       regionmeandat)
  zoomsnpdata<-rbind(zoomsnpdata,
                     regionsnpsdat)
}

colnames(zoompidat)<-c(names(zoompidat)[1:8],"region")
colnames(zoompimeandat)<-c(names(zoompimeandat)[1:4],"region")
colnames(zoomsnpdata)<-c(names(zoomsnpdata)[1:16],"region")
head(zoomsnpdata)

unique(zoompidat$region)
length(unique(zoompidat$region))/2

half<-length(unique(zoompimeandat$region))/2
end<-length(unique(zoompimeandat$region))

#Upper
ggplot()+
  geom_line(data=zoompidat[zoompidat$region %in% 
                             unique(zoompidat$region)[1:half],],
            aes(WinPos/1000000,pi,colour=rep),alpha=1/2)+
  geom_line(data=zoompimeandat[zoompimeandat$region %in% 
                                 unique(zoompidat$region)[1:half],],
                           aes(winpos/1000000,pi),size=0.5)+
  geom_point(data=zoomsnpdata[zoomsnpdata$region %in% 
                                unique(zoomsnpdata$region)[1:half],],
                              aes(Pos/1000000,0.01),
             shape=21,colour="black",fill="red")+
  geom_hline(yintercept = 0,colour="grey")+
  xlab("")+
  ylab(expression(pi))+
  scale_y_continuous(limits=c(0,0.035),breaks = c(0,0.015,0.03))+
  scale_colour_manual("Replicate",
                      values=pal)+
  facet_grid(tr~region,scales="free")+
  my.theme+
  theme(axis.text=element_text(size=10),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "top",
        legend.title = element_text(size=10,face="bold"),
        legend.text = element_text(size=10),
        legend.key = element_rect(fill="white",colour="grey"),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())
#Lower
ggplot()+
  geom_line(data=zoompidat[zoompidat$region %in% 
                             unique(zoompidat$region)[half:end],],
            aes(WinPos/1000000,pi,colour=rep),alpha=1/2)+
  geom_line(data=zoompimeandat[zoompimeandat$region %in% 
                                 unique(zoompidat$region)[half:end],],
            aes(winpos/1000000,pi),size=0.5)+
  geom_point(data=zoomsnpdata[zoomsnpdata$region %in% 
                                unique(zoomsnpdata$region)[half:end],],
             aes(Pos/1000000,0.01),
             shape=21,colour="black",fill="red")+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits=c(0,0.035),breaks = c(0,0.015,0.03))+
  scale_colour_manual("Replicate",
                      values=pal)+
  facet_grid(tr~region,scales="free")+
  my.theme+
  theme(axis.text=element_text(size=10),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "top",
        legend.title = element_text(size=10,face="bold"),
        legend.key = element_rect(fill="white",colour="grey"),
        legend.text = element_text(size=10),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())









ggplot()+
  geom_point(data=E_vs_M_pi,aes(Epi,Mpi,colour=Rep),alpha = 1/5)+
  xlab(expression(paste(pi," in E")))+
  ylab(expression(paste(pi," in M")))+
  facet_grid(Chr~Rep)+my.theme+
  theme(axis.text=element_text(size=10),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10))


