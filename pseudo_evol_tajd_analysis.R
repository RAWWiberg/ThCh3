##
# D.pseudoobscura project evolution linse
# PoPoolation tajd data analysis
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
E1D<-read.table("pseudo_evol_E1.D",header = FALSE,
                  sep = "\t")
colnames(E1D)<-c("Chr","WinPos","NrSNPs","PropGdReads","D")
E1D$D <- as.character(E1D$D)
# Subset to only known chromosomes
E1D<-E1D[grep("Unknown",E1D$Chr,invert = TRUE),]
E1D$tr<-rep("E",nrow(E1D))
E1D$rep<-rep("1",nrow(E1D))
str(E1D)
head(E1D)
nrow(E1D)
###
# M1
###
M1D<-read.table("pseudo_evol_M1.D",header = FALSE,
                 sep = "\t")
colnames(M1D)<-c("Chr","WinPos","NrSNPs","PropGdReads","D")
M1D$D <- as.character(M1D$D)
# Subset to only known chromosomes
M1D<-M1D[grep("Unknown",M1D$Chr,invert = TRUE),]
M1D$tr<-rep("M",nrow(M1D))
M1D$rep<-rep("1",nrow(M1D))
head(M1D)
nrow(M1D)
###
# E2
###
E2D<-read.table("pseudo_evol_E2.D",header = FALSE,
                 sep = "\t")
colnames(E2D)<-c("Chr","WinPos","NrSNPs","PropGdReads","D")
E2D$D <- as.character(E2D$D)
# Subset to only known chromosomes
E2D<-E2D[grep("Unknown",E2D$Chr,invert = TRUE),]
E2D$tr<-rep("E",nrow(E2D))
E2D$rep<-rep("2",nrow(E2D))
str(E2D)
head(E2D)
nrow(E2D)
###
# M2
###
M2D<-read.table("pseudo_evol_M2.D",header = FALSE,
                 sep = "\t")
colnames(M2D)<-c("Chr","WinPos","NrSNPs","PropGdReads","D")
M2D$D <- as.character(M2D$D)
# Subset to only known chromosomes
M2D<-M2D[grep("Unknown",M2D$Chr,invert = TRUE),]
M2D$tr<-rep("M",nrow(M2D))
M2D$rep<-rep("2",nrow(M2D))
head(M2D)
nrow(M2D)
###
# E3
###
E3D<-read.table("pseudo_evol_E3.D",header = FALSE,
                 sep = "\t")
colnames(E3D)<-c("Chr","WinPos","NrSNPs","PropGdReads","D")
E3D$D <- as.character(E3D$D)
# Subset to only known chromosomes
E3D<-E3D[grep("Unknown",E3D$Chr,invert = TRUE),]
E3D$tr<-rep("E",nrow(E3D))
E3D$rep<-rep("3",nrow(E3D))
str(E3D)
head(E3D)
nrow(E3D)
###
# M3
###
M3D<-read.table("pseudo_evol_M3.D",header = FALSE,
                 sep = "\t")
colnames(M3D)<-c("Chr","WinPos","NrSNPs","PropGdReads","D")
M3D$D <- as.character(M3D$D)
# Subset to only known chromosomes
M3D<-M3D[grep("Unknown",M3D$Chr,invert = TRUE),]
M3D$tr<-rep("M",nrow(M3D))
M3D$rep<-rep("3",nrow(M3D))
head(M3D)
nrow(M3D)
###
# E4
###
E4D<-read.table("pseudo_evol_E4.D",header = FALSE,
                 sep = "\t")
colnames(E4D)<-c("Chr","WinPos","NrSNPs","PropGdReads","D")
E4D$D <- as.character(E4D$D)
# Subset to only known chromosomes
E4D<-E4D[grep("Unknown",E4D$Chr,invert = TRUE),]
E4D$tr<-rep("E",nrow(E4D))
E4D$rep<-rep("4",nrow(E4D))
str(E4D)
head(E4D)
nrow(E4D)
###
# M4
###
M4D<-read.table("pseudo_evol_M4.D",header = FALSE,
                 sep = "\t")
colnames(M4D)<-c("Chr","WinPos","NrSNPs","PropGdReads","D")
M4D$D <- as.character(M4D$D)
# Subset to only known chromosomes
M4D<-M4D[grep("Unknown",M4D$Chr,invert = TRUE),]
M4D$tr<-rep("M",nrow(M4D))
M4D$rep<-rep("4",nrow(M4D))
head(M4D)
nrow(M4D)
# ####

# Combine D sets
# D for E and M as separate columns
# ####
E_vs_M_D<-cbind(as.character(c(M1D$D,M2D$D,M3D$D,M4D$D)),
                 as.character(c(E1D$D,E2D$D,E3D$D,E4D$D)),
                 c(rep("1",length(M1D$D)),
                   rep("2",length(M2D$D)),
                   rep("3",length(M3D$D)),
                   rep("4",length(M4D$D))),
                 c(as.character(E1D$Chr),
                   as.character(E2D$Chr),
                   as.character(E3D$Chr),
                   as.character(E4D$Chr)))

E_vs_M_D<-as.data.frame(E_vs_M_D)
colnames(E_vs_M_D)<-c("MD","ED","Rep","Chr")
#E_vs_M_D<-E_vs_M_D[E_vs_M_D$ED != "na" & E_vs_M_D$MD != "na",]
E_vs_M_D$MD<-as.numeric(as.character(E_vs_M_D$MD))
E_vs_M_D$ED<-as.numeric(as.character(E_vs_M_D$ED))
E_vs_M_D$Chr<-gsub("_group",".",E_vs_M_D$Chr)
head(E_vs_M_D)
str(E_vs_M_D)
nrow(E_vs_M_D)
# ####

# D for E an M (melted)
# ####
Ddat<-rbind(E1D,M1D,
             E2D,M2D,
             E3D,M3D,
             E4D,M4D)
# Remove "na"
#Ddat<-Ddat[Ddat$D != "na",]
Ddat$D<-as.numeric(Ddat$D)
Ddat$chrty<-vector(length=nrow(Ddat))
for (row in seq(1,nrow(Ddat))){
  if (grepl("X",as.character(Ddat$Chr[row]))){
    Ddat$chrty[row]<-"X"
    }
  else{
    Ddat$chrty[row]<-"A"
  }
}

# What is average diversity across all chromosomes
tapply(Ddat$D,INDEX = list(Ddat$tr,Ddat$chrty,Ddat$rep),
       mean,na.rm=TRUE)
# What is the error
tapply(Ddat$D,INDEX = list(Ddat$tr,Ddat$chrty,Ddat$rep),se)*1.69

# E lines seem to have lower diversity overall:
# indicative of selective sweeps (or BG selection)?

EDdat<-rbind(E1D,
              E2D,
              E3D,
              E4D)
#EDdat<-EDdat[EDdat$D != "na",]
EDdat$D<-as.numeric(EDdat$D)
EDdat$Chr<-factor(EDdat$Chr)
MDdat<-rbind(M1D,
              M2D,
              M3D,
              M4D)
#MDdat<-MDdat[MDdat$D != "na",]
MDdat$D<-as.numeric(MDdat$D)
MDdat$Chr<-factor(MDdat$Chr)
head(MDdat)
# Get mean D across treatments
###
# E
###
Emeandata<-tapply(EDdat$D,
                 INDEX = list(EDdat$WinPos,EDdat$Chr),mean,na.rm=TRUE)
head(Emeandata)
winpos<-dimnames(Emeandata)[[1]]
chrs<-dimnames(Emeandata)[[2]]
Emeandata<-as.data.frame(Emeandata)
Emeandata$WinPos<-rownames(Emeandata)
head(Emeandata)
str(Emeandata)

Ecombmeandata<-data.frame(winpos=vector(),D=vector,Chr=vector())
for(chr in chrs){
  D<-Emeandata[,chr]
  chr<-rep(chr,length(D))
  winpos<-as.numeric(rownames(Emeandata))
  data<-data.frame(winpos=winpos,D=D,Chr=chr)
  Ecombmeandata<-rbind(Ecombmeandata,data)
}
Ecombmeandata<-na.omit(Ecombmeandata)
Ecombmeandata$D<-as.numeric(Ecombmeandata$D)
rownames(Ecombmeandata)<-seq(1,nrow(Ecombmeandata))
Ecombmeandata$tr<-rep("E",nrow(Ecombmeandata))
str(Ecombmeandata)
head(Ecombmeandata)
unique(Ecombmeandata$Chr)

###
# M
###
Mmeandata<-tapply(MDdat$D,
                 INDEX = list(MDdat$WinPos,MDdat$Chr),mean,na.rm=TRUE)
str(Mmeandata)
winpos<-dimnames(Mmeandata)[[1]]
chrs<-dimnames(Mmeandata)[[2]]
Mmeandata<-as.data.frame(Mmeandata)
Mmeandata$WinPos<-rownames(Mmeandata)
head(Mmeandata)
str(Mmeandata)

Mcombmeandata<-data.frame(winpos=vector(),D=vector,Chr=vector())
for(chr in chrs){
  D<-Mmeandata[,chr]
  chr<-rep(chr,length(D))
  winpos<-as.numeric(rownames(Mmeandata))
  data<-data.frame(winpos=winpos,D=D,Chr=chr)
  Mcombmeandata<-rbind(Mcombmeandata,data)
}
Mcombmeandata<-na.omit(Mcombmeandata)
Mcombmeandata$D<-as.numeric(Mcombmeandata$D)
rownames(Mcombmeandata)<-seq(1,nrow(Mcombmeandata))
Mcombmeandata$tr<-rep("M",nrow(Mcombmeandata))
head(Mcombmeandata)

unique(Mcombmeandata$Chr)

Dmeandat<-rbind(Mcombmeandata,Ecombmeandata)
Dmeandat$Chr<-gsub("_group",".",Dmeandat$Chr)
head(Dmeandat)
head(Ddat)
Ddat$Chr<-gsub("_group",".",Ddat$Chr)
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
  geom_line(data=Ddat,
            aes(WinPos/1000000,D,colour=rep),alpha=1/2)+
  geom_line(data=Dmeandat,aes(winpos/1000000,D),size=0.5)+
  xlab("")+
  ylab(expression(D))+
  scale_colour_manual("Replicate",
                      values=pal)+
  geom_rect(data=chim_reg_dat,aes(xmin=start/1000000,
                                  xmax=end/1000000,ymin=0.03,ymax=0.035),
            fill = "black")+
  geom_text(data=chim_reg_dat,aes(x=(end+200000)/1000000,
                                  y=0.03,label=reg))+
  facet_grid(Chr~tr)+
  my.theme+
  theme(axis.text=element_text(size=10),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "top",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

# Extract the data for each chimney region
chim_reg_dat
head(Ddat)
head(Dmeandat)
zoomDdat<-data.frame()
zoomDmeandat<-data.frame()
zoomsnpdata<-data.frame()
for(reg in chim_reg_dat$reg){
  chr<-as.character(chim_reg_dat$Chr[chim_reg_dat$reg==reg])
  start<-chim_reg_dat$start[chim_reg_dat$reg==reg]-50000
  end<-chim_reg_dat$end[chim_reg_dat$reg==reg]+50000
  region<-paste(chr,":",start/1000000,"-",end/1000000,sep="")
  
  regiondat<-Ddat[Ddat$Chr==chr & Ddat$WinPos>start & 
                     Ddat$WinPos<end,]
  regiondat<-cbind(regiondat,
                   rep(region,nrow(regiondat)))
  
  regionmeandat<-Dmeandat[Dmeandat$Chr==chr & Dmeandat$winpos>start & 
                             Dmeandat$winpos<end,]
  regionmeandat<-cbind(regionmeandat,
                       rep(region,nrow(regionmeandat)))
  regionsnpsdat<-glmdat_top_q[glmdat_top_q$Chr2==chr & 
                                glmdat_top_q$Pos > start & 
                                glmdat_top_q$Pos < end,]
  regionsnpsdat<-cbind(regionsnpsdat,
                       rep(region,nrow(regionsnpsdat)))
  zoomDdat<-rbind(zoomDdat,
                   regiondat)
  zoomDmeandat<-rbind(zoomDmeandat,
                       regionmeandat)
  zoomsnpdata<-rbind(zoomsnpdata,
                     regionsnpsdat)
}

colnames(zoomDdat)<-c(names(zoomDdat)[1:8],"region")
colnames(zoomDmeandat)<-c(names(zoomDmeandat)[1:4],"region")
colnames(zoomsnpdata)<-c(names(zoomsnpdata)[1:16],"region")

unique(zoomDdat$region)
length(unique(zoomDdat$region))/2

half<-length(unique(zoomDmeandat$region))/2
end<-length(unique(zoomDmeandat$region))

#Upper
ggplot()+
  geom_hline(yintercept = 0,colour="grey")+
  geom_line(data=zoomDdat[zoomDdat$region %in% 
                             unique(zoomDdat$region)[1:half],],
            aes(WinPos/1000000,D,colour=rep),alpha=1/2)+
  geom_line(data=zoomDmeandat[zoomDmeandat$region %in% 
                                 unique(zoomDdat$region)[1:half],],
                           aes(winpos/1000000,D),size=0.5)+
  geom_point(data=zoomsnpdata[zoomsnpdata$region %in% 
                                unique(zoomsnpdata$region)[1:half],],
                              aes(Pos/1000000,2),
             shape=21,colour="black",fill="red")+
  scale_y_continuous(limits=c(-2,2),breaks = c(-2,-1,0,1,2))+
  xlab("")+
  ylab(expression(paste("Tajima's ",italic(D))))+
  scale_colour_manual("Replicate",
                      values=pal)+
  facet_grid(tr~region,scales="free")+
  my.theme+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11,face="bold"),
        strip.text.x = element_text(angle=0,size=12),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "top",
        legend.title = element_text(size=10,face="bold"),
        legend.key = element_rect(fill="white",colour="grey"),
        legend.text = element_text(size=10),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())


#Lower
ggplot()+
  geom_line(data=zoomDdat[zoomDdat$region %in% 
                             unique(zoomDdat$region)[half+1:end],],
            aes(WinPos/1000000,D,colour=rep),alpha=1/2)+
  geom_line(data=zoomDmeandat[zoomDmeandat$region %in% 
                                 unique(zoomDdat$region)[half+1:end],],
            aes(winpos/1000000,D),size=0.5)+
  geom_point(data=zoomsnpdata[zoomsnpdata$region %in% 
                                unique(zoomsnpdata$region)[half+1:end],],
             aes(Pos/1000000,2),
             shape=21,colour="black",fill="red")+
  geom_hline(yintercept = 0,colour="grey")+
  scale_y_continuous(limits=c(-2,2),breaks = c(-2,-1,0,1,2))+
  xlab("")+
  ylab("")+
  scale_colour_manual("Replicate",
                      values=pal,guide=FALSE)+
  facet_grid(tr~region,scales="free")+
  my.theme+
  theme(axis.text=element_text(size=10),
        strip.text.x = element_text(angle=0,size=12),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "top",
        legend.title = element_text(size=10,face="bold"),
        legend.key = element_rect(fill="white",colour="grey"),
        legend.text = element_text(size=10),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())









ggplot()+
  geom_point(data=E_vs_M_D,aes(ED,MD,colour=Rep),alpha = 1/5)+
  xlab(expression(paste(D," in E")))+
  ylab(expression(paste(D," in M")))+
  facet_grid(Chr~Rep)+my.theme+
  theme(axis.text=element_text(size=10),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10))


