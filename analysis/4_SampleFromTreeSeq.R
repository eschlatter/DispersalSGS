#setwd("/projectnb/dispevol/E_Schlatter/kernelSGS/analysis")
setwd("C:/Users/eschlatter/Dropbox/DispersalSGS/analysis")
treefile = "../output/500k/ts_8349801080707846925_t500"
library(tidyverse)
library(vcfR)

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
# read in empirical sample sites data
#############################################################################
sample_sites <- read_csv('../data/sgs_sites_SLiM',show_col_types=FALSE)
sample_inds <- data.frame(site=numeric(),id=numeric(),dist=numeric(),x=numeric(),y=numeric()) # to store the sampled individuals

# then run Fst, using a list of the lists of Node ID's as the sample_sets
for(site in 1:nrow(sample_sites)){
  if(sample_sites$n_SNP[site]>0){ # if empirical sample there is >0
    site_x <- sample_sites$x[site]
    site_y <- sample_sites$y[site]
    inds$dist <- sqrt((inds$x-site_x)^2+(inds$y-site_y)^2)
    #  inds[[paste0('site',site,'dist')]] <- sqrt((inds$x-site_x)^2+(inds$y-site_y)^2)
    inds_for_site <- slice_min(inds,order_by=dist,n=sample_sites$n_SNP[site])
    temp.df <- data.frame(site=sample_sites$Site[site],
                          id=inds_for_site$id,
                          dist=inds_for_site$dist,
                          x = inds_for_site$x,
                          y = inds_for_site$y)
    sample_inds <- rbind(sample_inds,temp.df)  
  }
}

which(duplicated(sample_inds$id)) # check for individuals included in more than one sample; there are none

# output:

# save this for later reference of which individuals were sampled
write.csv(sample_inds,file=paste0(treefile,"_sample_inds.csv")) 
# vector of individual IDs
individs=as.vector(sample_inds$id)
write.table(individs, file=paste0(treefile,"_individuals.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
# vector of individual ID names ("tsk_i")
sample_inds <- mutate(sample_inds,gt_name=paste0("tsk_",id))
write.table(sample_inds$gt_name,file=paste0(treefile,"_individualnames.txt"),
            eol=" ",quote=FALSE,row.names=FALSE,col.names=FALSE)