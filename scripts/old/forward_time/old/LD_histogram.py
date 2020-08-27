# set population size
# set sample size
# set unlinked loci
# eset mutation rate
# initialise population
# evolve for n generations
# check ld
# estimate size

import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from itertools import combinations
import numpy as np


Ne = 100
S = 50
n_loci = 100
gens = 3
repeats = 500
import matplotlib.pyplot as plt

r2s = []
estimates = []
for r in range(repeats):
    print(r+1)

    # set up population
    pop = sim.Population(size=[Ne], ploidy=2, loci=[1]*n_loci)

    # evolve population
    pop.evolve( 
        initOps = [sim.InitSex(), sim.InitGenotype(freq = [0.5,0.5])],
        matingScheme = sim.RandomMating(),
        gen = gens
    )

    # take sample of size S
    sample_pop = drawRandomSample(pop, sizes = S)

    # get allele frequency
    sim.stat(sample_pop, alleleFreq = range(0,n_loci), vars = ['alleleFreq'])
    allele_freqs = sample_pop.vars()['alleleFreq']
    seg_alleles = []
    # find which alleles are segregating
    for k in allele_freqs.keys():
        if (allele_freqs[k][0] > 0.04) and (allele_freqs[k][1] > 0.04):
            seg_alleles.append(k)

    # only proceed if there are 2 or more segregating alleles (to measure r2)
    if len(seg_alleles) < 2:
        pass

    # calculate r2 values
    sim.stat(sample_pop, LD = list(combinations(seg_alleles, r=2 )), vars = ['R2'])

    # calculate mean r2
    r2_total = 0
    count = 0
    for layer1 in sample_pop.vars()['R2'].values():
        for value in layer1.values():
            #print(value)
            r2_total += value
            count +=1
    mean_r2 = r2_total / count

    # correct r2 for sample size
    r2_drift = (mean_r2 - 1/(2*S)) / (1 - 1/(2*S))

    #get Ne estimate
    Ne_est = 1/(3*r2_drift)
    estimates.append(Ne_est)
    r2s.append(r2_drift)

print(np.mean(estimates))

plt.hist(estimates, range = [-Ne, 3*Ne], bins = 300)
plt.plot([Ne,Ne], [0,20])
plt.xlim(-Ne,3*Ne)
plt.show()

print(r2s)
plt.hist(r2s, bins = 10)
plt.plot([1/(3*Ne),1/(3*Ne)], [0,20])
plt.xlim(0, .02)
plt.show()