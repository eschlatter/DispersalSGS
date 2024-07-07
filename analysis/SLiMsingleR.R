#  Created by Sam Champer, 2020.
#  A product of the Messer Lab, http://messerlab.org/slim/

#  Sam Champer, Ben Haller and Philipp Messer, the authors of this code, hereby
#  place the code in this file into the public domain without restriction.
#  If you use this code, please credit SLiM-Extras and provide a link to
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.
#  Thank you.

# Running SLiM in R is easy and clean!

# Run slim using the system() function. Store the output in slim_out.
slim_out <- system("slim -d N_sponge=50000 -d F_goby=50 -d NMAX=0.25 -d MATECHOICE_DIST=0.1 2_SimulateEvolution.slim", intern=TRUE)
print(slim_out)