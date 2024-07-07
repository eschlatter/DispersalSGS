#  Created by Sam Champer, 2020.
#  A product of the Messer Lab, http://messerlab.org/slim/

#  Sam Champer, Ben Haller and Philipp Messer, the authors of this code, hereby
#  place the code in this file into the public domain without restriction.
#  If you use this code, please credit SLiM-Extras and provide a link to
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.
#  Thank you.

# Running SLiM in R is easy and clean!

# A little more advanced:
# Running an array of simulations in parallel:

library(foreach)
library(doParallel)
library(future)
library(parallelly)

cl <- makeCluster(future::availableCores())
registerDoParallel(cl)

# Vectors to vary the params
NSPONGE <- c(50000,50000,50000,50000)
FGOBY <- c(60,70,50,50)
N_MAX <- c(0.25,0.2,0.2,0.2)
MATEDIST <- c(0.1,0.1,0.2,0.15)

# Run SLiM in parallel:
raw_slim_output_matrix <- foreach(i=1:length(NSPONGE)) %dopar% {
    # Use string maniputaion functions to configure the command line args,
    # then run SLiM with system(),
    # then keep only the last line.
    slim_out <- system(sprintf("slim -d N_sponge=%f -d F_goby=%f -d NMAX=%f -d MATECHOICE_DIST=%f 2_SimulateEvolution.slim",
                               NSPONGE[i], FGOBY[i],N_MAX[i],MATEDIST[i]), intern=TRUE)
  }
stopCluster(cl)


# 
# Elements in raw_slim_output_matrix contain all of the text that SLiM outputs.
# It could be parsed and put into a matrix, or just flattened into a list, as below:
parsed_output = c()
for (col in raw_slim_output_matrix)
  for (row in col)
    for (line in row)
      if (startsWith(line, "OUT:"))
        # Strip off the "OUT:" from the beginning of the line.
        parsed_output = c(parsed_output, substring(line, first=5))

save(parsed_output,file='../output/cluster_output_test.RData')