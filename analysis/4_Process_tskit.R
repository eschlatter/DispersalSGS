#setwd("/projectnb/dispevol/E_Schlatter/kernelSGS/analysis")
setwd("C:/Users/eschlatter/Dropbox/DispersalSGS/analysis")
treefile = "../output/ts_8238198510834684392_t160000"

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
inds <- mutate(inds,x=str_trim(str_sub(x,start=18)),pedID=str_trim(str_sub(pedID,-8,-1)),
               pedp1=str_trim(str_sub(pedp1,-8,-1)),pedp2=str_trim(str_sub(pedp2,-8,-1)),
               age=str_trim(str_sub(age,-1,-1)),sex=str_trim(str_sub(sex,-1,-1)))
inds <- mutate(inds,id=0:(nrow(inds)-1),x=as.numeric(x),y=as.numeric(y),age=as.numeric(age))

#############################################################################
# read in sample sites data
#############################################################################
sample_sites <- read_csv('../data/sgs_sites_SLiM',show_col_types=FALSE)

# create dataframe to store the sampled individuals
sample_inds <- data.frame(site=numeric(),id=numeric(),dist=numeric(),x=numeric(),y=numeric())

# and fill it
for(site in 1:nrow(sample_sites)){
  if(sample_sites$n_SNP[site]>0){ # skip any sites with no sampling
    
    # sampling site location
    site_x <- sample_sites$x[site]
    site_y <- sample_sites$y[site]
    
    # calculate distance of each individual from the sampling site location
    inds$dist <- sqrt((inds$x-site_x)^2+(inds$y-site_y)^2)
    
    # pick the closest individuals
    inds_for_site <- slice_min(inds,order_by=dist,n=sample_sites$n_SNP[site])
    
    # store data about them, and add it to the running dataframe
    temp.df <- data.frame(site=sample_sites$Site[site],
                          id=inds_for_site$id,
                          dist=inds_for_site$dist,
                          x = inds_for_site$x,
                          y = inds_for_site$y)
    sample_inds <- rbind(sample_inds,temp.df)  
  }
}

if(length(which(duplicated(sample_inds$id)))>0){
  print(paste('Warning: same individuals sampled at multiple sites!', which(duplicated(sample_inds$id))))
}

# save a csv of the sampled individuals, for later reference
write.csv(sample_inds,file=paste0(treefile,"_sample_inds.csv"))

#############################################################################
# read in VCF file from tskit
#############################################################################

simVCF <- read.vcfR(paste0(treefile,".vcf"),verbose=FALSE)

# keep just the columns of @gt that match sample_inds$id
sample_inds <- mutate(sample_inds,gt_name=paste0("tsk_",id))
col_keep <- colnames(simVCF@gt) %in% sample_inds$gt_name
col_keep[1] <- TRUE #keep first column of @gt (called FORMAT)
simVCF@gt <- simVCF@gt[,col_keep]

# create a factor of the sampling site of each individual in vcf
ind_sites <- as.factor(arrange(sample_inds,id)$site)

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
ld_report <- gl.report.ld.map(simGL,maf=0.02)
simGL <- gl.filter.ld(simGL,ld_report, threshold=0.8) # R^2>0.8 is what was used in the paper
simGL@pop <- ind_sites #restore site info

# rare alleles
simGL <- gl.filter.maf(simGL,threshold=0.02) # also removes monomorphic loci

#############################################################################
# Calculate pairwise fsts
#############################################################################

fsts <- gl.fst.pop(simGL,nboots=1)
save(fsts,file=paste0(treefile,"_FST.RData"))