import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
from allelic import get_migration_matrix


Ne = 100
S = 50
n_loci = 50
gens = 20
repeats = 50
n_subpops = 10
ms = [0.001,0.01,0.05,0.1,0.2,.5]
initial_frequencies = [.5, .5]




###########################


def get_mean_r2(Ne, S, n_loci, gens, repeats, n_subpops,initial_frequencies, m):
    M = get_migration_matrix(m, n_subpops)
    pop = sim.Population(size=[Ne]*n_subpops, ploidy=2, loci=[1]*n_loci, infoFields = 'migrate_to')
    sim.initGenotype(pop, freq = initial_frequencies)
    pop.evolve( 
            initOps = [sim.InitSex(), sim.InitGenotype(freq = initial_frequencies)],
            preOps = sim.Migrator(rate=M),
            matingScheme = sim.RandomMating(),
            gen = gens
        )
    
    sample_pop = drawRandomSample(pop, sizes = [S]*n_subpops)

    # get allele frequencies
    sim.stat(sample_pop, alleleFreq = range(0,n_loci), vars = ['alleleFreq_sp'])
    # calculate r2 values
    sim.stat(sample_pop, LD = list(combinations(list(range(n_loci)), r=2)), vars = ['R2_sp'])

    r2s = []

    for sp in range(n_subpops):

        allele_freqs = sample_pop.dvars(sp).alleleFreq
        seg_alleles = [k for k in range(n_loci) if np.abs(.5-allele_freqs[k][0]) < .5 - 0.05]
        if len(seg_alleles) < 2: raise Exception("<2 segregating alleles")

        r2_sum = count = 0

        for pairs in combinations(seg_alleles, r = 2):
            r2_sum += sample_pop.dvars(sp).R2[pairs[0]][pairs[1]]
            count += 1
        
        mean_r2 = r2_sum / count
        r2s.append(mean_r2)

    return r2s


def get_r2_drift(r2_tot, S):
    return np.array(r2_tot) - 1/(2*S) / (1 - 1/(2*S))


def get_Ne(r2_drift):
    return 1/(3*r2_drift)



Ne_m = {}
for m in ms:
    print(m)
    Ne_m[m] = []
    for r in range(repeats):
        print(r+1)
        r2_tot = get_mean_r2(Ne, S, n_loci, gens, repeats, n_subpops,initial_frequencies, m)
        r2_drift = get_r2_drift(r2_tot, S)
        Ne_ests = get_Ne(r2_drift)
        Ne_m[m].append(Ne_ests)

print(Ne_m)

###########################

means = [np.mean(Ne_m[m]) for m in ms]
plt.scatter(ms, means, edgecolors='black', color = 'white')
plt.plot([min(ms),max(ms)], [100,100], 'k--')

# plot housekeeping
plt.plot([min(ms), max(ms)], [Ne,Ne], 'k--')
plt.ylim(Ne/2,max(means)*1.1)
plt.xlim(.9*min(ms), 1.1*max(ms))
plt.xscale("log")
plt.xticks(ticks = ms,labels = ms)
plt.xlabel('Migration rate ')
plt.ylabel(r'$\hat{N}_e$')


plt.show()



# print(len(estimates) / np.sum(1.0 / np.array(estimates)))
# plt.hist(estimates, range = [-Ne, 3*Ne], bins = 300)
# plt.plot([Ne,Ne], [0,20])
# plt.xlim(-Ne,3*Ne)
# plt.show()

# print(r2s)
# plt.hist(r2s, bins = 10)
# plt.plot([1/(3*Ne),1/(3*Ne)], [0,20])
# plt.xlim(0, .02)
# plt.show()