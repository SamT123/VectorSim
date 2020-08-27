
import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt

Ne = 100
S = 50
n_loci = 1000
gens = 3
repeats = 10
n_subpops = 10
ms = [0.01,0.05,0.1,0.2,.5, .9]
initial_frequencies = [0.9,0.1]


##########################

def get_FCs():
    ''''Runs simulations for allelic fluctuations model with n subpopulations, and returns a list of FC values (one for each subpopulation)'''



def fc_variant(freqs_init, freqs_end):
    xi_1, yi_1 = freqs_init[0], freqs_end[0]
    xi_2, yi_2 = freqs_init[1], freqs_end[1]

    return (1/2) * (xi_1-yi_1)**2 / ( (xi_1 + yi_1)/2 - xi_1*yi_1)  + (1/2) * (xi_2-yi_2)**2 / ( (xi_2 + yi_2)/2 - xi_2*yi_2)




##########################

all_fcs = {}
all_Nes = {}
for m in ms:
    m_adj = m / (n_subpops-1)
    M = np.full( (n_subpops,n_subpops), m_adj )
    np.fill_diagonal(M, 0)
    M = M.tolist()

    Ne_ests = []
    FCs = []
    for r in range(repeats):
        print(r+1)



        popNe = sim.Population(size=[Ne]*n_subpops, ploidy=2, loci=[1]*n_loci, infoFields = 'migrate_to')
        popS  = sim.Population(size=[S], ploidy=2, loci=[1]*n_loci)
        sim.initGenotype(popNe, freq = initial_frequencies)
        sim.initGenotype(popS , freq = initial_frequencies)

        sim.stat(popS, alleleFreq = range(0,n_loci), vars = ['alleleFreq'])
        initial_allele_freqs = popS.vars()['alleleFreq']


        popNe.evolve( 
            initOps = [sim.InitSex()],
            preOps = sim.Migrator(rate=M),
            matingScheme = sim.RandomMating(),
            gen = gens
        )

        sample_pop = drawRandomSample(popNe, sizes = [S]+[0]*9)
        sim.stat(sample_pop, alleleFreq = range(0,n_loci), vars = ['alleleFreq'])


        meanFC_repeat = 0
        countFC = 0


        allele_freqs = sample_pop.vars()['alleleFreq']

        for k in allele_freqs.keys() :
            inits = list(initial_allele_freqs[k].values())
            ends = list(allele_freqs[k].values())

            if len(inits) + len(ends) == 4:
                meanFC_repeat += fc_variant(list(initial_allele_freqs[k].values()), list(allele_freqs[k].values()))
                countFC += 1

        if countFC != 0:
            meanFC_repeat = meanFC_repeat / countFC

            Ne_repeat = gens / (2*(meanFC_repeat - 1/(2*S) - 1/(2*S) ))

            Ne_ests.append(Ne_repeat)
            FCs.append(meanFC_repeat)

    all_fcs[m] = FCs
    all_Nes[m] = Ne_ests





FCmean_Nes = [gens / (2*(np.mean(all_fcs[m]) - 1/(2*S) - 1/(2*S) )) for m in ms] 
print(FCmean_Nes)
plt.plot(ms, FCmean_Nes)
plt.show()
