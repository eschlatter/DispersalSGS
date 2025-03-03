#setwd("/projectnb/dispevol/E_Schlatter/kernelSGS/analysis")
setwd("C:/Users/eschlatter/Dropbox/DispersalSGS/analysis")
#treefile <- commandArgs(TRUE)[1]
#treefile = "C:/Users/eschlatter/Dropbox/DispersalSGS/output/checkdir/1_1639670251/ts_1639670251_t3"
treefile = "C:/Users/eschlatter/Dropbox/DispersalSGS/output/test_fst_metrics_50k/1_2138423022/ts_2138423022_t10"
treefile

library(dplyr)
library(stringr)
library(readr)
library(tidyr)

library(vcfR)
library(dartR)
library(snpStats)


#############################################################################
# read in Individuals data from tskit
#############################################################################
inds <- read_csv(paste0(treefile,".csv"),col_names=F,show_col_types=FALSE)
inds <- as.data.frame(t(inds))
inds <- separate_wider_delim(inds,cols=V1,delim=',',names=c('flags','x','y','zero','parents','array','dtype','pedID','pedp1','pedp2','age','subpop','sex','flags2')) %>%
  dplyr::select(x,y,pedID,pedp1,pedp2,age,sex)

inds <- mutate(inds,x = as.numeric(gsub("[^0-9.]","",x))) %>%
  mutate(y = as.numeric(gsub("[^0-9.]","",y))) %>%
  mutate(pedID = as.factor(gsub("[^0-9.]","",pedID))) %>%
  mutate(pedp1 = str_sub(gsub("\\D","",pedp1),start=2)) %>%
  mutate(pedp2 = str_sub(gsub("\\D","",pedp2),start=2)) %>%
  mutate(age = as.numeric(gsub("\\D","",age))) %>%
  mutate(sex = as.numeric(gsub("\\D","",sex))) %>%
  mutate(id=0:(nrow(inds)-1))

#############################################################################
# read in Sites data from tskit
#############################################################################

sites <- read_csv(paste0(treefile,"_sites.csv"),col_names=F)
sites <- as.data.frame(t(sites))
colnames(sites) = c('pedID','site_member')
sites <- mutate(sites,pedID=as.factor(pedID))


#############################################################################
# read in VCF file from tskit
#############################################################################

simVCF <- read.vcfR(paste0(treefile,".vcf"),verbose=FALSE)


#############################################################################
# connect genetic info to site info
#############################################################################

# make a column in inds so that ids exactly match the vcf
inds <- mutate(inds,gt_name=paste0("tsk_",id))
inds <- left_join(inds,sites,by="pedID")
# create a factor of the sampling site of each individual in vcf
ind_sites <- as.factor(arrange(inds,id)$site_member)

#############################################################################
# Convert to genind, add population info, convert to genlight
#############################################################################

simGI <- vcfR2genind(simVCF)
simGI@pop <- ind_sites
simGL <- gi2gl(simGI,verbose=0)

#############################################################################
# Filter: LD and rare alleles
#############################################################################

# LD:
simGL@pop <- as.factor(rep('pop1',length(simGL@ind.names))) #remove site info temporarily
ld_report <- gl.report.ld.map(simGL,maf=0.02, plot.out=FALSE)
simGL <- gl.filter.ld(simGL,ld_report, threshold=0.8, verbose=3) # R^2>0.8 is what was used in the paper
simGL@pop <- ind_sites #restore site info

# rare alleles
simGL <- gl.filter.maf(simGL,threshold=0.02, plot.out=FALSE, verbose=3) # also removes monomorphic loci

#############################################################################
# Calculate pairwise fsts
#############################################################################

fsts <- gl.fst.pop(simGL,nboots=1, verbose=3)
save(fsts,file=paste0(treefile,"_FST.RData"))