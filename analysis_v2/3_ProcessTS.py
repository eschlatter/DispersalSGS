import msprime, tskit, pyslim, pandas, csv, sys, numpy

#####################################################
# read in tree sequence from SLiM and add genetics
#####################################################

# treefile = sys.argv[1]
treefile = "C:/Users/eschlatter/Dropbox/DispersalSGS/output/checkdir/1_1641294773/ts_1641294773_t3"
ts = tskit.load("%s.trees" % (treefile))
max_num_roots = max([t.num_roots for t in ts.trees()])
print(f"Maximum number of roots: {max_num_roots}")

# recapitate
ts_r = pyslim.recapitate(ts, ancestral_Ne=420000, recombination_rate=1.5e-8)

# add mutations
ts_m = msprime.sim_mutations(
    ts_r,
    model=msprime.SLiMMutationModel(type=0, next_id=pyslim.next_slim_mutation_id(ts_r)),
    rate=5e-9
)

# convert to nucleotides
nts = pyslim.convert_alleles(pyslim.generate_nucleotides(ts_m)) # randomly generate nucleotides

#####################################################
# select only sampled individuals to save
#####################################################

t = nts.tables
slim_meta = nts.metadata["SLiM"]["user_metadata"]
v_pedID = slim_meta["ids"]
v_sites = slim_meta["tags"]
subset_inds = numpy.nonzero(v_sites)[0]
v_pedID_sub = [v_pedID[i] for i in subset_inds]
v_sites_sub = [v_sites[i] for i in subset_inds]

# pull out the individuals whose pedigree ID is in the subset of v_pedID indexed by subset_inds
ind_table_indices = []
for indiv_i in range(0,len(t.individuals)-1):
    indiv = t.individuals[indiv_i]
    if indiv.metadata["pedigree_id"] in v_pedID_sub:
        ind_table_indices.append(indiv_i)

# pull out the ids of nodes with individual in ind_table_indices
node_table_indices = []
for node_i in range(0, len(t.nodes)-1):
    if t.nodes[node_i].individual in ind_table_indices:
        node_table_indices.append(node_i)

# pass those nodes to simplify to get a subsetted tree sequence
subset_ts = nts.simplify(node_table_indices)

#####################################################
# output
#####################################################

# 1) VCF file (by tskit ID)
with open(("%s.vcf" % (treefile)),"w") as f:
    subset_ts.write_vcf(f)

# 2) CSV with tskit ID and pedigree ID, so we can link them
    # (also contains location info, sex, etc)
with open(("%s.csv" % (treefile)),"wt",newline='') as testfile:
    csv_writer = csv.writer(testfile, delimiter=',')
    csv_writer.writerow(subset_ts.tables.individuals)

# 3) CSV with pedigree ID and sampling site
site_info = [v_pedID_sub, v_sites_sub]
with open(("%s_sites.csv" % (treefile)),"wt",newline='') as testfile:
    csv_writer = csv.writer(testfile, delimiter=',')
    csv_writer.writerows(site_info)
