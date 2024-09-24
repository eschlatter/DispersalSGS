#setwd("/projectnb/dispevol/E_Schlatter/kernelSGS/analysis")
setwd("C:/Users/eschlatter/Dropbox/DispersalSGS/analysis")
treefile = "../output/500k/ts_8349801080707846925_t500"
library(tidyverse)
library(vcfR)

# read in locations of simulated individuals
sample_inds <- read_csv(paste0(treefile,"_sample_inds.csv")) %>%
  mutate(site=as.factor(site))

# get VCF file containing the sample of simulated individuals
simVCF <- read.vcfR(paste0(treefile,".vcf"),verbose=FALSE)

# convert to genind, and add populations
simGI <- vcfR2genind(simVCF) #vcfR
simGI@pop=as.factor(sample_inds$site)
#print(simGI)

# convert to genlight
simGL <- gi2gl(simGI,verbose=0) #dartR

#############################################################################
# Filter: LD and minor alleles
#############################################################################

# LD:
simGL@pop <- as.factor(rep('pop1',length(simGL@ind.names))) #remove site info temporarily
ld_report <- gl.report.ld.map(simGL,maf=0.02)
simGL <- gl.filter.ld(simGL,ld_report, threshold=0.8) # R^2>0.8 is what was used in the paper. But they used PLINK, which does a sliding window thing...(?)
simGL@pop <- as.factor(sample_inds$site) #restore site info

# HWE
# I think I'm going to skip this. It's normally meant to remove loci under selection, and there is no selection in this simulation. Probably what it's actually doing is removing loci with signs of population subdivision, which I specifically don't want to do.
# simGL <- gl.filter.hwe(simGL,subset='each')

# rare alleles
simGL <- gl.filter.maf(simGL,threshold=0.02) # also removes monomorphic loci
# after filtering for rare alleles: 1118 genotypes, 219 binary SNPs

#############################################################################
# Calculate pairwise fsts
#############################################################################

fst_sim <- gl.fst.pop(simGL,nboots=1)  #dartR
fst_sim <- as.data.frame(fst_sim) %>%
  rownames_to_column(var = 'site1')
fst_sim <- pivot_longer(fst_sim,cols = 2:ncol(fst_sim),
                        names_to = c('site2'))
fst_sim$sites <- NA
for (j in 1:nrow(fst_sim)) {
  fst_sim$sites[j] <- paste(sort(as.matrix(fst_sim[j, 1:2])), collapse = '.')
}
fst_sim <- filter(fst_sim, !is.na(value)) %>%
  rename(sim_fst = value)
fst_sim <- arrange(fst_sim,sites)   
save(fst_sim,file=paste0(treefile,"_FST.RData"))
