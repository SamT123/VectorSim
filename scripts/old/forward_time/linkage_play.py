from linkage_functions import get_migration_matrix, get_mean_r2, get_mean_r2, get_Ne, get_r2_drift, hmean, calculate_estimates
import numpy as np
import matplotlib.pyplot as plt
import itertools
import simuPOP as sim
from simuPOP.sampling import drawRandomSample

def check_freq(a, crit):
    if a[0] > crit and a[1]>crit:
        return True
    return False
simu = False
if simu:
    m=.5
    S=200
    nLoci = 100
    pop = sim.Population(loci = [1]*nLoci, size=[500]*2, infoFields='migrate_to')
    pop.evolve(
        
        initOps=[sim.InitSex(),sim.InitGenotype(freq = [.5, .5])],
        preOps=sim.Migrator(rate=[[0,m],[m,.0]]),
        matingScheme=sim.RandomMating(),
        postOps=[
            sim.Stat(popSize=True),
            sim.PyEval('subPopSize'),
            sim.PyOutput('\n')
        ],

        gen = 100
    )
    sim.Migrator(rate=[[0,m],[m,.0]])

    sample_pop = drawRandomSample(pop, sizes = [S]*2)
    sim.stat(sample_pop, LD = list(itertools.combinations(list(range(nLoci)), r=2)), vars = ['R2_sp'])
    sim.stat(sample_pop, alleleFreq = range(0,nLoci), vars = ['alleleFreq_sp'])
    r2s = sample_pop.dvars().__dict__['subPop'][0]['R2']
    af = sample_pop.dvars().__dict__['subPop'][0]['alleleFreq']
    print(af)
    ks = list(r2s.keys())
    total = 0
    count = 0
    for L1 in ks:
        for L2 in list(r2s[L1].keys()):
            if check_freq(af[L1], .05) and check_freq(af[L2], .05):
                total += r2s[L1][L2]
                count += 1

    r2_mean = total / count
    r2_drift = r2_mean -  (1/(2*S) / (1 - 1/(2*S)))
    Ne = 1/(3*r2_drift)
    print(Ne)


############################
# now with msprime
############################

import msprime
import tskit

reps = 3
diploid_Ne = 100
diploid_S = 50

S = 2*diploid_S
mutation_rate = 5e-9 / 40
recom_rate = 1e-8
positions = [0, 1e8-1, 1e8, 2e8-1]
rates = [recom_rate, 0.5, recom_rate, 0]
num_loci = int(positions[-1])
n_subpops = 2
M = 1
m = M/n_subpops
population_configurations = [msprime.PopulationConfiguration(sample_size=S)] + [msprime.PopulationConfiguration(sample_size=0) for i in range(n_subpops - 1)]

migration_matrix = [
    [0,m],
    [m,0]]

recombination_map = msprime.RecombinationMap(
    positions=positions, rates=rates, num_loci=num_loci)
ests = []

for rep in range(reps):
    print(rep+1)
    tree_sequence = msprime.simulate(
        Ne=diploid_Ne, recombination_map=recombination_map,
        model="dtwf", mutation_rate = mutation_rate, 
        population_configurations=population_configurations,
        migration_matrix=migration_matrix)


    first_chrom2 = 'not assigned'
    n1=0
    n2=0
    for variant in tree_sequence.variants():
        if variant.position > 1e8:
            first_chrom2 = variant.index
            n2+=1
        if variant.position < (1e8 - 1):
            n1+=1

    r2 = tskit.LdCalculator(tree_sequence).r2_matrix()
    fil = np.zeros((n1+n2,n1+n2))
    fil[:n1, n1:] = 1

    i=0
    for v in tree_sequence.variants():
        if sum(v.genotypes / len(v.genotypes)) < 0.05:
            fil[:,i] = 0
            fil[i,:] = 0
        i += 1

    fil = fil.astype(int)
    r2 = np.where(fil, r2, np.zeros_like(r2))

    r2_mean = np.sum(r2) / np.sum(fil)
    r2_drift = r2_mean -  (1/(2*diploid_S) / (1 - 1/(2*diploid_S)))
    Ne_est = 1/(3*r2_drift)
    ests.append(Ne_est)
    print(n1+n2)
    print(Ne_est)

plt.hist(ests)
plt.plot([diploid_Ne,diploid_Ne], [0,5], color = "red")
plt.plot([n_subpops*diploid_Ne, n_subpops*diploid_Ne], [0,5], color = "blue")
plt.plot([np.mean(ests), np.mean(ests)], [0,5], color = "green")
plt.show()