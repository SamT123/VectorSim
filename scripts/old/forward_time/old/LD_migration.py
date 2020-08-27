import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt


Ne = 100
S = 50
n_loci = 50
gens = 200
repeats = 5
n_subpops = 10
ms = [0.01,0.05,0.1,0.2,.5, .9]





###########################
def get_mean_r2():
    

###########################

full_estimates = {}
for m in ms:
    m_adj = m / (n_subpops-1)
    M = np.full( (n_subpops,n_subpops), m_adj )
    np.fill_diagonal(M, 0)
    M = M.tolist()

    r2s = []
    estimates = []

    for r in range(repeats):
        print(r+1)

        # set up population
        pop = sim.Population(size=[Ne]*n_subpops, ploidy=2, loci=[1]*n_loci, infoFields = 'migrate_to')

        # evolve population
        pop.evolve( 
            initOps = [sim.InitSex(), sim.InitGenotype(freq = [0.5,0.5])],
            preOps = sim.Migrator(rate=M),
            matingScheme = sim.RandomMating(),
            gen = gens
        )

        # take sample of size S
        sample_pop = drawRandomSample(pop, sizes = [S]*n_subpops)

        # get allele frequency
        sim.stat(sample_pop, alleleFreq = range(0,n_loci), vars = ['alleleFreq_sp'])

        # calculate r2 values
        sim.stat(sample_pop, LD = list(combinations(list(range(n_loci)), r=2)), vars = ['R2_sp'])

        estimates.append([])
        r2s.append([])
        for sp in range(n_subpops):
            allele_freqs = sample_pop.dvars(sp).alleleFreq

            seg_alleles = []
            # find which alleles are segregating
            for k in allele_freqs.keys():
                if (allele_freqs[k][0] > 0.04) and (allele_freqs[k][1] > 0.04):
                    seg_alleles.append(k)

            # only proceed if there are 2 or more segregating alleles (to measure r2)
            if len(seg_alleles) < 2:
                continue



            # calculate mean r2
            r2_total = 0
            count = 0

            for pairs in combinations(seg_alleles, r=2):

                r2_i = sample_pop.dvars(sp).R2[pairs[0]][pairs[1]]
                r2_total += r2_i
                count+=1


            mean_r2 = r2_total / count

            # correct r2 for sample size
            r2_drift = (mean_r2 - 1/(2*S)) / (1 - 1/(2*S))

            #get Ne estimate
            Ne_est = 1/(3*r2_drift)
            estimates[-1].append(Ne_est)
            r2s[-1].append(r2_drift)
    full_estimates[m] = estimates


means = [np.mean(full_estimates[m]) for m in ms]
plt.scatter(ms, means, edgecolors='black', color = 'white')
plt.plot([min(ms),max(ms)], [100,100], 'k--')
plt.xscale('log')
plt.xticks(ticks = ms, labels = ms)
plt.ylim(50,150)
plt.xlim(min(ms)*0.95,max(ms)*1.05)
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