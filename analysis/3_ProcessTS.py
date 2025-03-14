import msprime, tskit, pyslim, pandas, csv

treefile = "../output/500k/ts_8349801080707846925_t500"
#treefile = "../output/ts_8238198510834684392_t160000"
ts = tskit.load("%s.trees" % (treefile))
max_num_roots = max([t.num_roots for t in ts.trees()])
print(f"Maximum number of roots: {max_num_roots}")

# recapitate
ts_r = pyslim.recapitate(ts, ancestral_Ne=42000, recombination_rate=1e-8)

# add mutations
ts_m = msprime.sim_mutations(
    ts_r,
    model=msprime.SLiMMutationModel(type=0, next_id=pyslim.next_slim_mutation_id(ts_r)),
    rate=5e-9
)
print(f"Number of mutations: {ts_m.num_mutations}")

# convert to nucleotides
nts = pyslim.convert_alleles(pyslim.generate_nucleotides(ts_m)) # randomly generate nucleotides

#Write to VCF
with open(("%s.vcf" % (treefile)),"w") as f:
    nts.write_vcf(f)
print(f"VCF written")

# write 'individuals' table to csv
t = nts.tables
indy = t.individuals
with open(("%s.csv" % (treefile)),"wt",newline='') as testfile:
    csv_writer = csv.writer(testfile, delimiter=',')
    csv_writer.writerow(indy)