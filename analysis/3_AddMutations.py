import msprime, tskit, pyslim

ts = tskit.load("../output/slim_2003563054.trees")

# individuals have a "location" property:
ind = ts.individual(0)
ind



#ts = ts.simplify()

# check for coalescence
for t in ts.trees():
    assert t.num_roots == 1, ("not coalesced! on segment {} to {}".
        format(t.interval[0], t.interval[1]))

#mutated = msprime.mutate(ts, rate=1e-7, random_seed=1, keep=True)
#mutated.dump("../output/final_overlaid.trees")
