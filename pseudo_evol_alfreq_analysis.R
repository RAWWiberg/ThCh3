##
# D.pseudoobscura project evolution linse
# PoPoolation2 allele frequencies data analysis
# Last Modified: March 2017
##
#clean environment
rm(list = ls(all = TRUE))
#
setwd("~/Desktop/Data/RData/")

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
source("RScripts/ggplot_theme.R")
options(scipen = 6)

# Read in data
alfreq<-read.table("pseudo_evol_trcomp_pwc_red_chr.tab",
                  sep = "\t")
#colnames(alfreq)<-c(names(alfreq)[1:8],"R1","R2","R3","R4")
colnames(alfreq)<-c("chrom","pos","ref","alleles","states",
                    "V6","type","maj","R1","R2","R3","R4")
str(alfreq)
names(alfreq)
head(alfreq)


# How many multiallelic
nrow(alfreq[alfreq$alleles > 2,])
nrow(alfreq[alfreq$alleles == 2,])
# Subset to biallelic
alfreq<-alfreq[alfreq$alleles == 2,]
# How many sites
nrow(alfreq)
# How many sites are fixed in all comparisons? (diff = 1)
nrow(alfreq[alfreq$R1 == 1 & 
              alfreq$R2 == 1 & 
              alfreq$R3 == 1 &
              alfreq$R4 == 1,])

# How many sites show at least one fixation
nrow(alfreq[alfreq$R1 == 1 | 
              alfreq$R2 == 1 | 
              alfreq$R3 == 1 |
              alfreq$R4 == 1,])

# Count the number of fixed differences in windows
# Load chromosome info
chrom_dat<-read.table("~/Desktop/Data/pseudo_project/reference/dpse_r31_FB2013_02/pseudo_chrom_data.tab",
                      header = TRUE, sep = ",")

head(chrom_dat[chrom_dat$Chr == "2",])
head(alfreq)

# Melt the data
alfreq_m<-melt(data = alfreq,id.vars = c("chrom","pos"),
               measure.vars = c("R1","R2","R3","R4"))
colnames(alfreq_m)<-c("chr","pos","rep","diff")
head(alfreq_m)

nrow(alfreq[alfreq$pos >= 1 & alfreq$pos <= 100000 & alfreq$R1 == 1,])
nrow(alfreq[alfreq$pos >= 1 & alfreq$pos <= 100000 & alfreq$R2 == 1,])
nrow(alfreq[alfreq$pos >= 1 & alfreq$pos <= 100000 & alfreq$R3 == 1,])
nrow(alfreq[alfreq$pos >= 1 & alfreq$pos <= 100000 & alfreq$R4 == 1,])

# Make some vectors
dxy_v<-vector()
nsubs_v<-vector()
npoly_v<-vector()
avg_diff_v<-vector()
wind_v<-vector()
chr_v<-vector()
start_v<-vector()
end_v<-vector()
nsnps_v<-vector()
rep_v<-vector()
win_len<-50000
win_ovlap<-0
for(rep in unique(alfreq_m$rep)){
  for(chr in unique(alfreq_m$chr)){
    max_l<-chrom_dat$length[as.character(chrom_dat$Chr) == chr]
    chr_dat<-alfreq_m[as.character(alfreq_m$chr) == chr &
                        as.character(alfreq_m$rep) == rep,]
    start <- 1
    end <- start + win_len
    wind <- 1
    while(end < max_l){
      chr_v <- c(chr_v,chr)
      start_v <- c(start_v,start)
      end_v <- c(end_v,end)
      wind_v <- c(wind_v, wind)
      nsnps<-nrow(chr_dat[chr_dat$pos >= start & 
                            chr_dat$pos <= end,])
      nsubs <- nrow(chr_dat[chr_dat$pos >= start & 
                            chr_dat$pos <= end & 
                            chr_dat$diff == 1,])
      npoly <- nrow(chr_dat[chr_dat$pos >= start & 
                            chr_dat$pos <= end &
                            chr_dat$diff != 1,])
      avg_diff<-mean(chr_dat$diff[chr_dat$pos >= start & 
                               chr_dat$pos <= end ])
      dxy <- nsubs/win_len
      dxy_v <- c(dxy_v,dxy)
      nsubs_v <- c(nsubs_v,nsubs)
      npoly_v <- c(npoly_v,npoly)
      avg_diff_v<-c(avg_diff_v,avg_diff)
      nsnps_v<-c(nsnps_v,nsnps)
      rep_v<-c(rep_v,rep)
      start<-end
      end <- start+win_len
      if(end > max_l){
        end <- max_l
      }
      wind<-wind+1
    }
  }
}
d_theta_dat<-data.frame("chr"=chr_v,
                     "start"=start_v,
                     "end"=end_v,
                     "wind"=wind_v,
                     "nsubs"=nsubs_v,
                     "dxy"=dxy_v,
                     "npoly"=npoly_v,
                     "avg_diff"=avg_diff_v,
                     "nsnps"=nsnps_v,
                     "rep"=rep_v)
d_theta_dat$WinPos <- (d_theta_dat$start+
                     ((d_theta_dat$end-d_theta_dat$start)/2))
d_theta_dat$chr<-gsub("_group",".",d_theta_dat$chr)
head(d_theta_dat)
summary(d_theta_dat)

pal <- c('#d7191c','#fdae61','#abd9e9','#2c7bb6')
# Plot
# Plot dXY
adiff_plot_dxy<-ggplot(data=d_theta_dat)+
  geom_line(aes(x = (start+((end-start)/2))/1000000,y=dxy,colour=rep))+
  #scale_y_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
  xlab("")+
  scale_x_continuous(breaks=c(0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30),
                     labels=c("0","","5","","10","","15","","20","","25","","30"))+
  ylab(expression(paste(italic(d)[XY])))+
  scale_colour_manual("Replicate",values=pal,labels=c("1","2","3","4"))+
  scale_y_continuous(breaks=c(0.000,0.003,0.006))+
  facet_grid(chr~.)

adiff_plot_dxy + my.theme + 
  theme(axis.text.x = element_text(size=10),
        axis.title = element_text(size=10),
        axis.text.y = element_text(size = 10,face = "bold"),
        strip.text.y = element_text(size = 10,angle = 0),
        legend.position = "top",
        legend.text = element_text(size=10,face="bold"),
        legend.title = element_text(size=10,face="bold"),
        legend.key = element_rect(fill="white",colour="grey"))

# Plot wat_theta
adiff_plot_wt<-ggplot(data=d_theta_dat)+
  geom_line(aes(x = (start+((end-start)/2)),y=wat_theta,colour=rep))+
  #scale_y_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
  xlab("Position")+
  ylab(expression(paste("Watterson's ",italic(theta))))+
  scale_colour_manual("Replicate",values=pal)+
  scale_x_continuous(breaks=c(0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30),
                     labels=c("0","","5","","10","","15","","20","","25","","30"))+
  scale_y_continuous(breaks=c(0.000,0.003,0.006))+
  facet_grid(chr~.)

adiff_plot_wt + my.theme + 
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size = 10,face = "bold"),
        strip.text.y = element_text(size = 10,angle = 0),
        axis.ticks.x = element_blank())





# Plot allele frequency differences, for "top SNPs"
# Get mean alfreq difference
meandata<-tapply(d_theta_dat$avg_diff,
                 INDEX = list(d_theta_dat$WinPos,d_theta_dat$chr),mean)

str(meandata)
winpos<-dimnames(meandata)[[1]]
chrs<-dimnames(meandata)[[2]]

combmeandata<-data.frame(winpos=vector(),adiff=vector,chr=vector())
for(chr in chrs){
  adiff<-meandata[,chr]
  chr<-rep(chr,length(adiff))
  winpos<-as.numeric(names(adiff))
  data<-data.frame(winpos=winpos,adiff=adiff,chr=chr)
  combmeandata<-rbind(combmeandata,data)
}
combmeandata<-na.omit(combmeandata)
rownames(combmeandata)<-seq(1,nrow(combmeandata))
str(combmeandata)
head(combmeandata)

adiff_plot<-ggplot(data=d_theta_dat)+
  geom_line(aes(x=(start+((end-start)/2))/1000000,
                 y=avg_diff,colour=rep),alpha=1/2)+
  geom_line(data=combmeandata,aes(x=winpos/1000000,
                y=adiff),colour="black",size=1)+
  scale_x_continuous(breaks=c(0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30),
                     labels=c("0","","5","","10","","15","","20","","25","","30"))+
  #scale_y_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
  xlab("")+
  ylab("Allele Frequency Differences")+
  scale_colour_manual("Replicate",values=pal)+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
  facet_grid(chr~.)

adiff_plot + my.theme + 
  theme(axis.text.x = element_text(size=10),
        axis.title = element_text(size=10),
        axis.text.y = element_text(size = 10,face = "bold"),
        strip.text.y = element_text(size = 10,angle = 0),
        legend.position = "top",
        legend.key = element_rect(fill="white",colour="grey"),
        legend.text = element_text(size=10,face="bold"),
        legend.title = element_text(size=10,face="bold"))


# Plot raw difference by SNP: histogram
head(alfreq_m)
rawdiff_plot<-ggplot(data=alfreq_m)+
  geom_histogram(aes(diff,fill=rep),position="dodge",alpha=1/2)+
  xlab("")+
  ylab("Allele Frequency Differences")+
  scale_fill_manual("Replicate",values=pal)

rawdiff_plot + my.theme + 
  theme(#axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12,face = "bold"),
    strip.text.y = element_text(size = 12,angle = 0),
    #axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size=10,face="bold"),
    legend.title = element_text(size=10,face="bold"))
