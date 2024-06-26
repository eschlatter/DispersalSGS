import msprime, tskit

ts = tskit.load("final.trees").simplify()
mutated = msprime.mutate(ts, rate=1e-7, random_seed=1, keep=True)
mutated.dump("final_overlaid.trees")