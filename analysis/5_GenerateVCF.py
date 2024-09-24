import msprime, tskit, pyslim, pandas, csv
treefile = "../output/500k/ts_8349801080707846925_t500"

with open(("%s_individuals.txt" % (treefile))) as f:
    individuals = f.readlines()
individuals = [int(numeric_string) for numeric_string in individuals]

with open(("%s_individualnames.txt" % (treefile))) as f:
    individualnames = f.readline()
individualnames=individualnames.split()
individualnames

nts = tskit.load("%s_temp.trees" % (treefile))

#Write to VCF
with open(("%s.vcf" % (treefile)),"w") as f:
    nts.write_vcf(f,individuals=individuals, individual_names=individualnames)