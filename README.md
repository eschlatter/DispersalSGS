# DispersalSGS

Using the neon goby Elacatinus lori as a model system, I am investigating whether contemporary dispersal patterns can be used to predict the spatial genetic structure (SGS) that has evolved over time. I will predict SGS by simulating the evolution of genomes moving through the contemporary seascape according to E. lori’s validated dispersal patterns, using the genomic simulation software SLiM. Then I will compare my simulated predictions to empirical SGS data to gain insights about the ecological processes that are most significant to the formation of SGS in a real, complex seascape. Finally, I will extend these methods to generate testable predictions about SGS in portions of E. lori’s range where those data are not available.

## Model Parameters and Assumptions
https://docs.google.com/spreadsheets/d/1zxMT583BmHr4OSTqSfwo9JgcTCxLh0YJXrepr66iRqI/edit?gid=0#gid=0

## Contents:
- `data`: folder containing raw data for reef mapping and locations of sampling sites

- `analysis`: folder containing scripts to perform data processing, simulation, and analysis
  - `1_GenerateMapAndParams.Rmd`: Preliminary data processing to generate map and parameter values for SLiM
  - `2_SimulateEvolution.slim`: SLiM simulation model
  - `3_AddMutations.py`: Process simulation output -- in particular, add neutral mutations.

- `output`: folder containing output files from analyses
    - `final.trees`: output of step 2: tree sequence recording of simulation
    - `final_overlaid.trees`: output of step 3: tree sequence with neutral mutations overlaid

- `map_reference.png`: A map of the study area and sampling sites, for visual reference
