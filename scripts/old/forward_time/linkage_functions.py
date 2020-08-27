import itertools
import numpy as np
import simuPOP as sim
import itertools
from simuPOP.sampling import drawRandomSample

def calculate_estimates(param_dict,  m, repeats, vary = False,):
    if vary:
        print('\t',vary[0], " simulations")
        param_dict[vary[0]] = vary[1]
        param_combos = list(itertools.product(param_dict['Ne'], param_dict['S'], param_dict['n_loci'], param_dict['t'], param_dict['n_subpops'], param_dict['initial_frequencies']))
        r2s = {}
        ests = {}
        ests_pool = {}
        for i in range(len(vary[1])):
            print(i)
            #print('\t\t = ', vary[1][i])
            r2s[vary[1][i]] = []
            ests[vary[1][i]] = []
            ests_pool[vary[1][i]] = []
            for _ in range(repeats):
                combo = param_combos[i]
                #print(combo)
                this_r2s = get_mean_r2(combo[0], combo[1], combo[2], combo[3], combo[4], combo[5], m)
                r2s[vary[1][i]].append(this_r2s)
                ests[vary[1][i]].append(get_Ne(get_r2_drift(np.array(this_r2s),param_dict['S'][0])) / combo[0])
            ests_pool[vary[1][i]] = get_Ne(get_r2_drift( np.mean(r2s[vary[1][i]] ),param_dict['S'][0])) / combo[0]

        return ests_pool,ests

    else:
        r2s = []
        ests = []
        ests_pool = []
        for i in range(repeats):
            print(i)
            #print(param_dict)
            this_r2s = get_mean_r2(param_dict['Ne'][0], param_dict['S'][0], param_dict['n_loci'][0], param_dict['t'][0], param_dict['n_subpops'][0], param_dict['initial_frequencies'][0], m)
            r2s.append(this_r2s)
            ests.append(get_Ne(get_r2_drift(np.array(this_r2s),param_dict['S'][0])) / param_dict['Ne'][0])
        ests_pool = get_Ne(get_r2_drift( np.mean(r2s),param_dict['S'][0])) / param_dict['Ne'][0]
        return ests_pool, ests



def hmean(arr):
    arr = np.array(arr).flatten()
    for i in arr:
        assert i > 0, "Negative estiamte!"

    return np.product(arr)**(1/len(arr))



def get_migration_matrix(m,n):

    if n == 1:
        return np.array([0])
    m_adj = m / (n)
    M = np.full( (n,n), m_adj )
    #np.fill_diagonal(M, 0)
    M = M.tolist()
    print(M)
    return M

def get_mean_r2(Ne, S, n_loci, gens, n_subpops,initial_frequencies, m):
    """Returns the mean r2 value for each subpopulation, in list of length n_subpops"""

    # make pairwise migration matrix
    M = get_migration_matrix(m, n_subpops)

    # initialise population
    n_alleles = 2
    pop = sim.Population(size=[Ne]*n_subpops, ploidy=2, loci=[1]*n_loci, alleleNames = [str(i) for i in range(n_alleles)], infoFields = 'migrate_to')
    sim.initGenotype(pop, freq = [initial_frequencies, 1-initial_frequencies])
    #sim.initGenotype(pop, freq = [1/n_alleles for i in range(n_alleles)])
    sim.initSex(pop)
    print(M)
    # run burn in generations
    pop.evolve( 
            initOps = [],
            preOps = sim.Migrator(M, mode = sim.BY_PROBABILITY),
            matingScheme = sim.RandomMating(),
            gen = gens
        )


    # take sample from each subpopulation
    sample_pop = drawRandomSample(pop, sizes = [S]+[0]*(n_subpops-1))
    #sim.dump(sample_pop)

    # get allele frequencies
    sim.stat(sample_pop, alleleFreq = range(0,n_loci), vars = ['alleleFreq_sp'])
    #print(sample_pop.dvars(0).alleleFreq)
    # calculate r2 values
    sim.stat(sample_pop, LD = list(itertools.combinations(list(range(n_loci)), r=2)), vars = ['R2_sp'])
    #print(sample_pop.dvars(0).R2)
    r2s = []

    for sp in [0]: #range(n_subpops*0):

        allele_freqs = sample_pop.dvars(sp).alleleFreq
        seg_alleles = [k for k in range(n_loci) if np.abs(.5-allele_freqs[k][0]) < .5 - 0.05]
        if len(seg_alleles) < 2: raise Exception("<2 segregating alleles")

        r2_sum = 0
        count = 0

        for pairs in itertools.combinations(seg_alleles, r = 2):
            r2_sum += sample_pop.dvars(sp).R2[pairs[0]][pairs[1]]
            count += 1
        
        mean_r2 = r2_sum / count
        r2s.append(mean_r2)

    return r2s


def get_r2_drift(r2_tot, S):
    return np.array(r2_tot) - 1/(2*S) / (1 - 1/(2*S))


def get_Ne(r2_drift):
    return 1/(3*r2_drift)
