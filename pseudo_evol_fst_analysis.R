##
# D.pseudoobscura project evolution linse
# PoPoolation2 fst data analysis
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
setwd("~/Desktop/Data/RData/")
source("RScripts/ggplot_theme.R")
options(scipen = 6)

# Read in data
R1fst<-read.table("pseudo_evol_R1_chr.fst",header = FALSE,
                  sep = "\t")
colnames(R1fst)<-c("Chr","WinPos","NrSNPs","Prop>min-cov","AvgMinCov","PWfst")
R1fst$dat<-rep("R1",nrow(R1fst))
R1fst$PWfst<-gsub("1:2=","",R1fst$PWfst)
R1fst$PWfst<-as.numeric(R1fst$PWfst)
R1fst$Chr2 <- factor(gsub("_group",".",R1fst$Chr))
R1fst <- R1fst[R1fst$'Prop>min-cov' > 0.5,]
head(R1fst)
summary(R1fst)


R2fst<-read.table("pseudo_evol_R2_chr.fst",header = FALSE,
                  sep = "\t")
colnames(R2fst)<-c("Chr","WinPos","NrSNPs","Prop>min-cov","AvgMinCov","PWfst")
R2fst$dat<-rep("R2",nrow(R2fst))
R2fst$PWfst<-gsub("1:2=","",R2fst$PWfst)
R2fst$PWfst<-as.numeric(R2fst$PWfst)
R2fst$Chr2 <- factor(gsub("_group",".",R2fst$Chr))
R2fst <- R2fst[R2fst$'Prop>min-cov' > 0.5,]
head(R2fst)
summary(R2fst)

R3fst<-read.table("pseudo_evol_R3_chr.fst",header = FALSE,
                  sep = "\t")
colnames(R3fst)<-c("Chr","WinPos","NrSNPs","Prop>min-cov","AvgMinCov","PWfst")
R3fst$dat<-rep("R3",nrow(R3fst))
R3fst$PWfst<-gsub("1:2=","",R3fst$PWfst)
R3fst$PWfst<-as.numeric(R3fst$PWfst)
R3fst$Chr2 <- factor(gsub("_group",".",R3fst$Chr))
R3fst <- R3fst[R3fst$'Prop>min-cov' > 0.5,]
head(R3fst)
summary(R3fst)

R4fst<-read.table("pseudo_evol_R4_chr.fst",header = FALSE,
                  sep = "\t")
colnames(R4fst)<-c("Chr","WinPos","NrSNPs","Prop>min-cov","AvgMinCov","PWfst")
R4fst$dat<-rep("R4",nrow(R4fst))
R4fst$PWfst<-gsub("1:2=","",R4fst$PWfst)
R4fst$PWfst<-as.numeric(R4fst$PWfst)
R4fst$Chr2 <- factor(gsub("_group",".",R4fst$Chr))
R4fst <- R4fst[R4fst$'Prop>min-cov' > 0.5,]
head(R4fst)
nrow(R4fst)
summary(R4fst)

# Combine datasets
combdat<-rbind(R1fst,R2fst,R3fst,R4fst)
summary(combdat)
nrow(combdat)
head(combdat)

# Sort first by chromosome then by position
combdat <- sort_df(combdat, vars = c("Chr2","WinPos"))
combdat$Chr <- factor(gsub("_.*","",combdat$Chr)) 
tapply(combdat$PWfst,
       INDEX = list(combdat$Chr2,combdat$dat),median)
tapply(combdat$PWfst,
       INDEX = list(combdat$Chr2,combdat$dat),mean)
combdat$Chr3<-gsub("\\..*","",combdat$Chr2)
combdat$Chr3<-gsub("X.*","X",combdat$Chr3)
unique(combdat$Chr3)

# Get mean fst for each chromosome
fstdat<-tapply(combdat$PWfst,
                INDEX = list(combdat$dat,combdat$Chr3),mean)
reps<-dimnames(fstdat)[[1]]#Replicate
chrs<-dimnames(fstdat)[[2]]#Chromosome
fstdat<-data.frame(fstdat)
colnames(fstdat)<-as.character(chrs)
fstdat$reps<-reps
fstdat<-melt(fstdat,id.vars = c("reps"))
colnames(fstdat)<-c("reps","chromosome","fst")
fstdat$chrt<-vector(length=nrow(fstdat))
for(i in 1:nrow(fstdat)){
  if(fstdat$chromosome[i] == "X"){
    fstdat$chrt[i]<-"X"
  }else{
    fstdat$chrt[i]<-"Autosome"
  }
}
# Split the data by Autosomes and X chromosome
A_fstdat<-fstdat[fstdat$chrt=="Autosome",]
A_fstdat$chromosome<-factor(A_fstdat$chromosome)
X_fstdat<-fstdat[fstdat$chrt=="X",]
X_fstdat$chromosome<-factor(X_fstdat$chromosome)
A_fstdatm<-tapply(A_fstdat$fst,INDEX=list(A_fstdat$reps),mean)
A_fstdatsd<-tapply(A_fstdat$fst,INDEX=list(A_fstdat$reps),sd)
X_fstdatm<-tapply(X_fstdat$fst,INDEX=list(X_fstdat$reps),mean)
X_fstdatsd<-tapply(X_fstdat$fst,INDEX=list(X_fstdat$reps),sd)

# Function for expected FST on X from:
# Ramachandran et al., 2004 Human Genomics 1:87-97

# FX is given by the following equation:
Fx<-function(faut,z){
  1-((9*(z+1)*(1-faut))/(8*((2*z)+1)-(1-faut)*((7*z)-1)))
}

# Get bootstrapped expected FST on X
head(combdat)
boots<-1000
R1fst_A<-R1fst[grep("X",R1fst$Chr,invert=TRUE),]
R2fst_A<-R2fst[grep("X",R2fst$Chr,invert=TRUE),]
R3fst_A<-R3fst[grep("X",R3fst$Chr,invert=TRUE),]
R4fst_A<-R4fst[grep("X",R4fst$Chr,invert=TRUE),]
# z is the ratio of effective male population size to effective
# female population size
z<-1
exp_x<-vector()
exp_x_n<-vector()
for(boot in 1:boots){
  cat("Boot: ",boot,"\n")
  bootcombdat<-rbind(R1fst_A[sample(1:nrow(R1fst_A),
                                    size=nrow(R1fst_A),replace=TRUE),],
                     R2fst_A[sample(1:nrow(R2fst_A),
                                    size=nrow(R2fst_A),replace=TRUE),],
                     R3fst_A[sample(1:nrow(R3fst_A),
                                    size=nrow(R3fst_A),replace=TRUE),],
                     R4fst_A[sample(1:nrow(R4fst_A),
                                    size=nrow(R4fst_A),replace=TRUE),])
  bootcombdat$Chr3<-gsub("\\..*","",bootcombdat$Chr2)
  #Calculate expected fst
  bootcombdat$expfst<-Fx(bootcombdat$PWfst,z)
  exp_x_m<-tapply(bootcombdat$expfst,INDEX=list(bootcombdat$dat),mean)
  print(exp_x_m)
  exp_x<-c(exp_x,as.vector(exp_x_m))
  exp_x_n<-c(exp_x_n,dimnames(exp_x_m)[[1]])
}
head(exp_x,n=50)
exp_xdat<-data.frame(fst=exp_x,reps=exp_x_n)

exp_X_fstdatm<-tapply(exp_xdat$fst,INDEX=list(exp_xdat$reps),mean)
exp_X_fstdatsd<-2*tapply(exp_xdat$fst,INDEX=list(exp_xdat$reps),se)
ggplot()+geom_histogram(data=exp_xdat,aes(fst),binwidth = 0.005)+
  facet_grid(reps~.)+xlim(0,1)
# Get expected from raw data and combine
fstdat<-data.frame(fst=c(A_fstdatm,X_fstdatm),fstsd=c(A_fstdatsd,X_fstdatsd),
                   chrt=c(rep("Autosome",length(A_fstdatm)),
                          rep("Observed X",length(A_fstdatm))))
fstdat$reps<-c(dimnames(A_fstdatm)[[1]],dimnames(X_fstdatm)[[1]])

fstdat<-rbind(fstdat,data.frame(fst=Fx(A_fstdatm,z),
                                chrt=rep("Expected X",length=length(A_fstdatm)),
                                fstsd=as.vector(exp_X_fstdatsd),
                                reps=dimnames(A_fstdatm)[[1]]))
fstdat
if(z==1){
  pl<-"A)"
}else if(z==5){
  pl<-"B)"
}
# What is the expected FST on the X in each comparison
pal <- c('#d7191c','#fdae61','#abd9e9','#2c7bb6')
ggplot()+
  geom_bar(data=fstdat,aes(x=reps,y=fst,fill=chrt,group=chrt),
                  stat="identity",width=0.4,position="dodge")+
  geom_errorbar(data=fstdat,aes(x=reps,ymin=(fst-fstsd),
                                ymax=(fst+fstsd),group=chrt),
                width=0.1,position=position_dodge(0.36),colour="darkred")+
  scale_fill_manual("",values=c("grey80","grey70","grey60"))+
  scale_x_discrete(labels=c("1","2","3","4"))+
  xlab("Replicate")+ 
  ylim(0,1)+
  ggtitle(paste(pl," z = ",z,sep=""))+
  ylab(expression(italic(F)[ST]))+
  my.theme + 
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text.y = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12,angle = 0),
        plot.title = element_text(size=12,face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10,face="bold"),
        legend.title = element_blank())



# Get the average pw fst for each window across samples.
meandata<-tapply(combdat$PWfst,
                 INDEX = list(combdat$WinPos,combdat$Chr2),mean)
winpos<-dimnames(meandata)[[1]]
chrs<-dimnames(meandata)[[2]]
head(meandata)
meandata<-as.data.frame(meandata)
meandata$WinPos<-rownames(meandata)
head(meandata)
str(meandata)
combmeandata<-data.frame(winpos=vector(),fst=vector,Chr2=vector())
for(chr in chrs){
  fst<-meandata[,chr]
  chr<-rep(chr,length(fst))
  winpos<-as.numeric(rownames(meandata))
  data<-data.frame(winpos=winpos,fst=fst,Chr2=chr)
  combmeandata<-rbind(combmeandata,data)
}
combmeandata<-na.omit(combmeandata)
combmeandata$fst<-as.numeric(combmeandata$fst)
rownames(combmeandata)<-seq(1,nrow(combmeandata))
str(combmeandata)
head(combmeandata)
fstdat<-tapply(combmeandata$fst,
               INDEX = list(combmeandata$Chr2),mean)


# Mean/Median on autosomes
hist(combdat$PWf[which(substr(combdat$Chr2,1,1) != "X")])
mean(combdat$PWf[which(substr(combdat$Chr2,1,1) != "X")])
se(combdat$PWf[which(substr(combdat$Chr2,1,1) != "X")])

#Mean/Median on X chromosome
hist(combdat$PWf[which(substr(combdat$Chr2,1,1) == "X")])
mean(combdat$PWf[which(substr(combdat$Chr2,1,1) == "X")])
se(combdat$PWf[which(substr(combdat$Chr2,1,1) == "X")])

###
# Barplot of FST throughout autosomes and X for each pairwise comparison.
###


###
# PLOT FST THROUGHOUT THE GENOME
###
pal <- c('#d7191c','#fdae61','#abd9e9','#2c7bb6')
# Plotting
fst_plot <- ggplot()+
  geom_line(data=combdat,aes(WinPos/1000000,PWfst,colour = dat),
            alpha=1/2)+
  geom_line(data=combmeandata,aes(winpos/1000000,fst),size=1)+
  ylab(expression(paste("Pairwise ", italic(F)[ST])))+
  xlab("")+
  scale_x_continuous(breaks=c(0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30),
                     labels=c("0","","5","","10","","15","","20","","25","","30"))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  scale_colour_manual("Replicate",values = pal,labels=c("1","2","3","4"))+
  facet_grid(Chr2~.)
fst_plot + my.theme + 
  theme(axis.text.x = element_text(size=10),
        axis.title = element_text(size=10),
        axis.text.y = element_text(size = 10,face = "bold"),
        strip.text.y = element_text(size = 10,angle = 0),
        legend.position = "top",
        legend.text = element_text(size=10,face="bold"),
        legend.title = element_text(size=10,face="bold"),
        legend.key = element_rect(fill="white",colour="grey"))

###
# PLOT FST WITHIN CHIMNEYS
###
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
chim_reg_dat
# Load SNP data
glmdat_top_q <- read.table(
  "pseudo_evol/pseudo_evol_qbglm-unp_minc16_cov17-49_top_qval_snps.rout",
  sep = "", header = TRUE)
head(glmdat_top_q)


head(combdat)
head(combmeandata)
zoomfstdat<-data.frame()
zoomfstmeandat<-data.frame()
zoomsnpdata<-data.frame()
for(reg in chim_reg_dat$reg){
  chr<-as.character(chim_reg_dat$Chr[chim_reg_dat$reg==reg])
  start<-chim_reg_dat$start[chim_reg_dat$reg==reg]-50000
  end<-chim_reg_dat$end[chim_reg_dat$reg==reg]+50000
  region<-paste(chr,":",start/1000000,"-",end/1000000,sep="")
  
  regiondat<-combdat[combdat$Chr2==chr & combdat$WinPos>start & 
                       combdat$WinPos<end,]
  regiondat<-cbind(regiondat,
                   rep(region,nrow(regiondat)))
  
  regionmeandat<-combmeandata[
    combmeandata$Chr2==chr & combmeandata$winpos>start & 
                            combmeandata$winpos<end,]
  regionmeandat<-cbind(regionmeandat,
                       rep(region,nrow(regionmeandat)))
  regionsnpsdat<-glmdat_top_q[glmdat_top_q$Chr2==chr & 
                                glmdat_top_q$Pos > start & 
                                glmdat_top_q$Pos < end,]
  regionsnpsdat<-cbind(regionsnpsdat,
                       rep(region,nrow(regionsnpsdat)))
  zoomfstdat<-rbind(zoomfstdat,
                  regiondat)
  zoomfstmeandat<-rbind(zoomfstmeandat,
                      regionmeandat)
  zoomsnpdata<-rbind(zoomsnpdata,
                     regionsnpsdat)
}

colnames(zoomfstdat)<-c(names(zoomfstdat)[1:8],"region")
colnames(zoomfstmeandat)<-c(names(zoomfstmeandat)[1:3],"region")
colnames(zoomsnpdata)<-c(names(zoomsnpdata)[1:17],"region")

unique(zoomfstdat$region)
length(unique(zoomDdat$region))/2

half<-length(unique(zoomfstmeandat$region))/2
end<-length(unique(zoomfstmeandat$region))
head(zoomfstdat)
# Get dataset of mean fst
means<-unlist(tapply(combmeandata$fst,
                     INDEX = list(combmeandata$Chr2),mean))
means_chr2<-names(means)
means<-vector(length=length(unique(zoomfstdat$region)))
regions<-unique(zoomfstdat$region)

for(reg in unique(regions)){
  chr<-unlist(strsplit(reg,split = ":"))[1]
  means[which(regions==reg)]<-mean(combmeandata$fst[combmeandata$Chr2 == chr])
}
meansdat<-data.frame(means=means,region=regions)
meansdat
head(zoomfstmeandat)
head(zoomsnpdata)
ggplot()+
  geom_line(data=zoomfstdat[zoomfstdat$region %in% 
                            unique(zoomfstdat$region),],
            aes(WinPos/1000000,PWfst,colour=dat),alpha=1/2)+
  geom_line(data=zoomfstmeandat[zoomfstmeandat$region %in% 
                                unique(zoomfstdat$region),],
            aes(winpos/1000000,fst),size=0.5)+
  geom_point(data=zoomsnpdata[zoomsnpdata$region %in% 
                                unique(zoomsnpdata$region),],
             aes(Pos/1000000,1),
             shape=21,colour="black",fill="red")+
  geom_hline(yintercept = 0,colour="grey")+
  geom_hline(data=meansdat,aes(yintercept = means),
             colour="grey60",linetype="dashed")+
  scale_y_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
  xlab("")+
  ylab(expression(italic(F[ST])))+
  scale_colour_manual("Replicate",
                      values=pal)+
  facet_wrap(~region,scales="free")+
  my.theme+
  theme(axis.text=element_text(size=10),
        strip.text.x = element_text(angle=0,size=10),
        strip.text.y = element_text(angle=0,size=10),
        legend.position = "top",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())




##
# All Pairwise Analysis
##
# Load data
allpw_fst<-read.table("pseudo_evol_allpw_chr.fst",header = FALSE,sep = "\t")
head(allpw_fst)
head(data.frame(do.call('rbind', strsplit(
  as.character(allpw_fst[,6]),"[=,]"))))
#Rename columns
for(col in seq(6,length(colnames(allpw_fst)))){
  #split the column
  split<-data.frame(do.call('rbind', strsplit(
    as.character(allpw_fst[,col]),"[=,]")))
  allpw_fst[,col]<-as.numeric(as.character(split[,2]))
  names<-split[,1]
  if(length(unique(names))==1){
    colnames(allpw_fst)[col]<-as.character(unique(names))
  }
}
colnames(allpw_fst)<-c("Chr","WinPos","NrSNPs","Prop>min-cov","AvgMinCov",
                       colnames(allpw_fst)[6:length(colnames(allpw_fst))])
allpw_fst$Chr<-gsub("_group",".",allpw_fst$Chr)
head(allpw_fst)

# Split the data into all E pw comparisons and all M pairwise comparisons
Epops<-as.character(c(2,4,6,8))
Ecomps<-c("2:4","2:6","2:8","4:6","4:8","6:8")
E_allpw_fst<-allpw_fst[,c(1:6,grep(paste(Ecomps,collapse="|"),colnames(allpw_fst)))]
Ez<-5
Mpops<-as.character(c(1,3,5,7))
Mcomps<-c("1:3","1:5","1:7","3:5","3:7","5:7")
M_allpw_fst<-allpw_fst[,c(1:6,grep(paste(Mcomps,collapse="|"),colnames(allpw_fst)))]

head(E_allpw_fst)
head(M_allpw_fst)

# GET THE MEANS ACROSS COMPARISONS FOR EACH WINDOW ON XSOME
treats<-c("E","M")
aggregate_means<-vector()
aggregate_pop<-vector()
aggregate_chr<-vector()
aggregate_pos<-vector()
aggregate_tr<-vector()
for(tr in treats){
  if(tr == "E"){
    xsomes<-E_allpw_fst[grep("X",E_allpw_fst$Chr),]
    pops<-Epops
  }else if(tr == "M"){
    xsomes<-M_allpw_fst[grep("X",M_allpw_fst$Chr),]
    pops<-Mpops
  }
  for(snp in 1:nrow(xsomes)){
    for(pop in pops){
      #    cat("SNP: ",snp,"-",as.character(autosomes[snp,1]),
      #        "-",autosomes[snp,2],"Population: ",pop,"\n")
      popfsts<-vector()
      for(col in colnames(xsomes)){
        if(pop %in% unlist(strsplit(col,":"))){
          popfsts<-c(popfsts,xsomes[snp,which(colnames(xsomes)==col)])
        }
      }
      aggregate_means<-c(aggregate_means,mean(popfsts,na.rm=TRUE))
      aggregate_pop<-c(aggregate_pop,pop)
      aggregate_chr<-c(aggregate_chr,as.character(xsomes[snp,1]))
      aggregate_pos<-c(aggregate_pos,xsomes[snp,2])
      aggregate_pos<-c(aggregate_pos,xsomes[snp,2])
      aggregate_tr<-c(aggregate_tr,tr)
    }
  }
}
x_meandat<-data.frame(mfsts=aggregate_means,
                      pop=aggregate_pop,
                      chr=aggregate_chr,
                      pos=aggregate_pos,
                      tr=aggregate_tr)
head(x_meandat)
unique(x_meandat$tr)
# GET OVERALL MEAN FOR EACH CHROMOSOME (only one in the case of X) for each sample
# Remove the groups from the chromosome name to get averages across the whol
x_meandat$chr<-gsub(".\\..*","",x_meandat$chr)

pop_m_fst<-vector(
  length=length(unique(x_meandat$pop))*length(unique(x_meandat$chr)))
pop_m_chr<-vector(
  length=length(unique(x_meandat$pop))*length(unique(x_meandat$chr)))
pop_m_tr<-vector(
  length=length(unique(x_meandat$pop))*length(unique(x_meandat$chr)))
pops<-unique(x_meandat$pop)
i<-1
for(chr in unique(x_meandat$chr)){
  for(pop in pops){
    #cat("Chromosome: ",chr,"Population: ",pop," Element: ", i,"\n")
    x_meandat_sub<-x_meandat[x_meandat$chr==chr &
                             x_meandat$pop == pop,]
    tr<-as.character(unique(x_meandat_sub$tr))
    fsts<-vector()
    chrs<-vector()
    pop_m_fst[i]<-mean(x_meandat_sub$mfsts,na.rm=TRUE)
    pop_m_chr[i]<-chr
    pop_m_tr[i]<-tr
    i <-i + 1
  }
}
x_popsnp_m_fst<-data.frame(pop=rep(pops,length(unique(x_meandat$chr))),
                           mfst=pop_m_fst,
                           chr=pop_m_chr,
                           tr=pop_m_tr,
                           chrt=rep("Oserved: X",length(pops)))
x_popsnp_m_fst

# GET THE MEANS FOR AUTOSOMES
treats<-c("E","M")
aggregate_means<-vector()
aggregate_pop<-vector()
aggregate_chr<-vector()
aggregate_pos<-vector()
aggregate_tr<-vector()
for(tr in treats){
  if(tr == "E"){
    asomes<-E_allpw_fst[grep("X",E_allpw_fst$Chr,invert=TRUE),]
    pops<-Epops
  }else if(tr == "M"){
    asomes<-M_allpw_fst[grep("X",M_allpw_fst$Chr,invert=TRUE),]
    pops<-Mpops
  }
  for(snp in 1:nrow(asomes)){
    for(pop in pops){
      #    cat("SNP: ",snp,"-",as.character(autosomes[snp,1]),
      #        "-",autosomes[snp,2],"Population: ",pop,"\n")
      popfsts<-vector()
      for(col in colnames(asomes)){
        if(pop %in% unlist(strsplit(col,":"))){
          popfsts<-c(popfsts,asomes[snp,which(colnames(asomes)==col)])
        }
      }
      aggregate_means<-c(aggregate_means,mean(popfsts,na.rm=TRUE))
      aggregate_pop<-c(aggregate_pop,pop)
      aggregate_chr<-c(aggregate_chr,as.character(asomes[snp,1]))
      aggregate_pos<-c(aggregate_pos,asomes[snp,2])
      aggregate_pos<-c(aggregate_pos,asomes[snp,2])
      aggregate_tr<-c(aggregate_tr,tr)
    }
  }
}
a_meandat<-data.frame(mfsts=aggregate_means,
                      pop=aggregate_pop,
                      chr=aggregate_chr,
                      pos=aggregate_pos,
                      tr=aggregate_tr)
head(a_meandat)
unique(a_meandat$tr)
# GET OVERALL MEAN FOR EACH CHROMOSOME FOR EACH SAMPLE
# Remove the groups from the chromosome name to get averages across the whol
a_meandat$chr<-gsub("\\..*","",a_meandat$chr)
pop_m_fst<-vector(
  length=length(unique(a_meandat$pop))*length(unique(a_meandat$chr)))
pop_m_chr<-vector(
  length=length(unique(a_meandat$pop))*length(unique(a_meandat$chr)))
pop_m_tr<-vector(
  length=length(unique(a_meandat$pop))*length(unique(a_meandat$chr)))
pops<-unique(a_meandat$pop)
i<-1
for(chr in unique(a_meandat$chr)){
  for(pop in pops){
    #cat("Chromosome: ",chr,"Population: ",pop," Element: ", i,"\n")
    a_meandat_sub<-a_meandat[a_meandat$chr==chr &
                               a_meandat$pop == pop,]
    tr<-as.character(unique(a_meandat_sub$tr))
    fsts<-vector()
    chrs<-vector()
    pop_m_fst[i]<-mean(a_meandat_sub$mfsts,na.rm=TRUE)
    pop_m_chr[i]<-chr
    pop_m_tr[i]<-tr
    i <-i + 1
  }
}
a_popsnp_m_fst<-data.frame(pop=rep(pops,length(unique(x_meandat$chr))),
                           mfst=pop_m_fst,
                           chr=pop_m_chr,
                           tr=pop_m_tr,
                           chrt=rep("Autosomes",length(pops)))
a_popsnp_m_fst
# Summarise (get mean and SD across autosomes for each pop)
#Autosomes
a_mfst<-tapply(a_popsnp_m_fst$mfst,INDEX=list(a_popsnp_m_fst$pop),mean)
pops<-dimnames(a_mfst)[[1]]
a_sdfst<-tapply(a_popsnp_m_fst$mfst,INDEX=list(a_popsnp_m_fst$pop),sd)
#Xsomes
x_mfst<-tapply(x_popsnp_m_fst$mfst,INDEX=list(x_popsnp_m_fst$pop),mean)
dimnames(x_mfst)[[1]]
x_sdfst<-tapply(x_popsnp_m_fst$mfst,INDEX=list(x_popsnp_m_fst$pop),sd)

fst_dat<-data.frame(fst=c(a_mfst,x_mfst),sd=c(a_sdfst,x_sdfst),pop=seq(1,8),
                    tr=rep(c("M","E"),4),
                    chrt=rep(c("Autosome","Observed: X"),each=8),
                    rep=rep(c("1","2","3","4"),each=2))

fst_dat<-rbind(fst_dat,data.frame(fst=Fx(a_mfst,rep(c(1,Ez),4)),
                                  sd=rep(NA,8),pop=seq(1,8),
                                  tr=rep(c("M","E"),4),
                                  chrt=rep("Expected: X",8),
                                  rep=rep(c("1","2","3","4"),each=2)))
fst_dat

# DO SOME BOOTSTRAPPING
boot_exp_fstx<-vector()
boot_exp_pop<-vector()
boot_exp_tr<-vector()
i<-1
boots<-1000
for(tr in treats){
  if(tr == "E"){
    #Take only autosomes
    asomes<-E_allpw_fst[grep("X",E_allpw_fst$Chr,invert=TRUE),]
    pops<-Epops
    z<-Ez
    }
  else if(tr == "M"){
    #Take only autosomes
    asomes<-M_allpw_fst[grep("X",M_allpw_fst$Chr,invert=TRUE),]
    pops<-Mpops
    z<-1
  }
  # Run bootsampling
  for(boot in 1:boots){
    aggregate_means<-vector()
    aggregate_pop<-vector()
    aggregate_chr<-vector()
    aggregate_pos<-vector()
    aggregate_tr<-vector()
    asomes_sub<-asomes[sample(1:nrow(asomes),nrow(asomes),replace=TRUE),]
    asomes_sub$Chr<-gsub("\\..*","",asomes_sub$Chr)
    cat("Boot: ",boot,"TR: ",tr,"z: ",z,"\n")
    for(win in 1:nrow(asomes_sub)){
      for(pop in pops){
        popfsts<-vector()
        for(col in colnames(asomes_sub)){
          if(pop %in% unlist(strsplit(col,":"))){
            popfsts<-c(popfsts,asomes_sub[win,which(colnames(asomes_sub)==col)])
          }
        }
        aggregate_means<-c(aggregate_means,mean(popfsts,na.rm=TRUE))
        aggregate_pop<-c(aggregate_pop,pop)
        aggregate_chr<-c(aggregate_chr,as.character(asomes_sub[win,1]))
        aggregate_pos<-c(aggregate_pos,asomes_sub[win,2])
        aggregate_pos<-c(aggregate_pos,asomes_sub[win,2])
        aggregate_tr<-c(aggregate_tr,tr)
      }
    }
    a_meandat<-data.frame(mfsts=aggregate_means,
                          pop=aggregate_pop,
                          chr=aggregate_chr,
                          pos=aggregate_pos,
                          tr=aggregate_tr)
    pop_fsts<-tapply(a_meandat$mfsts,INDEX=list(a_meandat$pop),mean,na.rm=TRUE)
    dimnames(pop_fsts)[[1]]
    
    print(pop_fsts)
    boot_exp_fstx<-c(boot_exp_fstx,Fx(as.vector(pop_fsts),z))
    boot_exp_pop<-c(boot_exp_pop,dimnames(pop_fsts)[[1]])
    boot_exp_tr<-c(boot_exp_tr,rep(tr,length(as.vector(pop_fsts))))
    i<-i+1
  }
}
bootsdata<-data.frame(exp_fstx=boot_exp_fstx,
                      pop=boot_exp_pop,
                      tr=boot_exp_tr)
head(bootsdata)
nrow(bootsdata)
ggplot()+
  geom_histogram(data=bootsdata,aes(exp_fstx),binwidth=0.005)+
  xlim(0,1)+
  facet_grid(pop~.)
# Get 2*se for the bootstrap data and add to fst_data
expfstx_sd<-2*tapply(bootsdata$exp_fstx,INDEX=list(bootsdata$pop),se)
dimnames(expfstx_sd)[[1]]
fst_dat$sd[fst_dat$chrt=="Expected: X"]<-as.vector(expfstx_sd)
fst_dat
# PLOT
pal <- c('#d7191c','#fdae61','#abd9e9','#2c7bb6')
ggplot()+
  geom_bar(data=fst_dat,aes(rep,fst,fill=chrt),
           position="dodge",stat="identity")+
  geom_errorbar(data=fst_dat,aes(rep,ymin=fst-sd,ymax=fst+sd,group=chrt),
                width=0.1,position=position_dodge(0.9),colour="darkred")+
  xlab("Replicate")+
  ylab(expression(italic(F)[ST]))+
  scale_fill_manual("",values=c("grey80","grey70","grey60"))+
  facet_grid(tr~.)+
  my.theme + 
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text.y = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12,angle = 0),
        plot.title = element_text(size=12,face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10,face="bold"),
        legend.title = element_blank())

# What is the mean autosome fst
cat("M")
mean(fst_dat$fst[
    fst_dat$chrt=="Autosome" & fst_dat$tr == "M"])
2*se(fst_dat$fst[
  fst_dat$chrt=="Autosome" & fst_dat$tr == "M"])

cat("E")
mean(fst_dat$fst[
  fst_dat$chrt=="Autosome" & fst_dat$tr == "E"])
2*se(fst_dat$fst[
  fst_dat$chrt=="Autosome" & fst_dat$tr == "E"])

# What is the ratio in E vs M
cat("M")
mean(fst_dat$fst[
  fst_dat$chrt=="Observed: X" & fst_dat$tr == "M"]/fst_dat$fst[
    fst_dat$chrt=="Autosome" & fst_dat$tr == "M"])
2*se(fst_dat$fst[
  fst_dat$chrt=="Observed: X" & fst_dat$tr == "M"]/fst_dat$fst[
    fst_dat$chrt=="Autosome" & fst_dat$tr == "M"])
cat("E")
mean(fst_dat$fst[
  fst_dat$chrt=="Observed: X" & fst_dat$tr == "E"]/fst_dat$fst[
    fst_dat$chrt=="Autosome" & fst_dat$tr == "E"])
2*se(fst_dat$fst[
  fst_dat$chrt=="Observed: X" & fst_dat$tr == "E"]/fst_dat$fst[
    fst_dat$chrt=="Autosome" & fst_dat$tr == "E"])










A<-seq(0.25,0.8,0.01)
eX_z1<-Fx(A,1)
eX_z2<-Fx(A,2)
eX_z3<-Fx(A,3)
eX_z4<-Fx(A,4)
eX_z5<-Fx(A,5)
expFST<-data.frame(Afst=rep(A,5),
                   eX=c(eX_z1,eX_z2,eX_z3,eX_z4,eX_z5),
                   z=rep(c("1","2","3","4","5"),each=length(A)))

ggplot()+
  geom_line(data=expFST,aes(Afst,eX,colour=z),size=1)+
  xlab(expression(paste("Observed Autosomal ",italic(F)[ST],sep="")))+
  ylab(expression(paste("Expected X chromosome ",italic(F)[ST],sep="")))+
  scale_color_grey()+ my.theme + 
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(size=12),
        axis.text.y = element_text(size = 12,face = "bold"),
        strip.text.y = element_text(size = 12,angle = 0),
        plot.title = element_text(size=12,face="bold"),
        legend.position = "right",
        legend.text = element_text(size=10,face="bold"))


plot(A,eX_z1,type="l")
points(A,eX_z5,type="l",col="blue")
points(A,eX_z2,type="l",col="green")
points(A,eX_z2,type="l",col="red")


