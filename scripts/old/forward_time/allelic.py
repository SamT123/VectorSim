import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import sys

##########################
def get_migration_matrix(m,n):

    if n == 1:
        return np.array([0])
    m_adj = m / (n)
    M = np.full( (n,n), m_adj )
    #np.fill_diagonal(M, 0)
    M = M.tolist()
    return M

def repair(frequencies):
    for allele in [0,1]:
        if allele not in frequencies.keys():
            frequencies[allele] = 0
    return frequencies

def get_FCs(Ne, S, n_loci, gens, n_subpops, initial_frequencies, m):
    ''''Runs simulations for allelic fluctuations model with n subpopulations, and returns a list of FC values (one for each subpopulation)'''
    # population to evolve ((from infinite gamete pool))
    popNe = sim.Population(size=[Ne]*n_subpops, ploidy=2, loci=[1]*n_loci, infoFields = 'migrate_to')
    # initial sample population (from infinite gamete pool)
    popS  = sim.Population(size=[S]*n_subpops, ploidy=2, loci=[1]*n_loci)
    sim.initGenotype(popNe, freq = initial_frequencies)
    sim.initGenotype(popS , freq = initial_frequencies)

    # get initial sample allele frequencies
    sim.stat(popS, alleleFreq = range(n_loci), vars = ['alleleFreq_sp'])

    M = get_migration_matrix(m, n_subpops)

    popNe.evolve( 
        initOps = [sim.InitSex()],
        preOps = sim.Migrator(rate=M),
        matingScheme = sim.RandomMating(),
        gen = gens
    )

    sample_pop = drawRandomSample(popNe, sizes = [S]*n_subpops)
    sim.stat(sample_pop, alleleFreq = range(n_loci), vars = ['alleleFreq_sp'])
    all_FCs = []

    for sp in range(n_subpops):
        initial_allele_frequencies = popS.dvars(sp).alleleFreq
        final_allele_frequencies = sample_pop.dvars(sp).alleleFreq
        sp_count = 0
        sp_FC = 0
        for locus in range(n_loci):
            init_pair = repair(initial_allele_frequencies[locus])
            end_pair = repair(final_allele_frequencies[locus])
            if init_pair[0]**2+init_pair[1]**2 != 1:
                sp_FC += fc_variant( [init_pair[0], init_pair[1]], [end_pair[0], end_pair[1]])
                sp_count += 1
        
        all_FCs.append(sp_FC / sp_count)

    return all_FCs


def fc_variant(freqs_init, freqs_end):
    xi_1, yi_1 = freqs_init[0], freqs_end[0]
    xi_2, yi_2 = freqs_init[1], freqs_end[1]

    if ( (xi_1 + yi_1)/2 - xi_1*yi_1) == 0:
        print(xi_1, yi_1)
    
    if ( (xi_2 + yi_2)/2 - xi_2*yi_2) == 0:
        print(xi_2, yi_2)

    return (1/2) * (xi_1-yi_1)**2 / ( (xi_1 + yi_1)/2 - xi_1*yi_1)  + (1/2) * (xi_2-yi_2)**2 / ( (xi_2 + yi_2)/2 - xi_2*yi_2)


def get_ests(FC_value, S, t):
    return t / (2* (FC_value - 1/(2*S) - 1/(2*S)) ) 

##########################
def main():
    Ne = 100
    Ne_list = [75,100,150, 200]
    S = 50
    S_list = [20,40,60]
    n_loci = 200
    loc_list = [50, 100, 200, 400]
    t = 3
    t_list = [3,5,7,9]
    repeats = 200
    n_subpops = 20
    n_subpops_list = [10,20,50]
    one_subpop = 1
    ms = [0.001, 0.01,0.05,0.1,0.5]
    no_m = 0
    initial_frequencies = [0.5,0.5]
    freq_list = [[0.5,0.5],[0.75,0.25], [0.95,0.05]]


    # is there bias without migration?

    ## Varying Ne
    print("Ne")
    Ne_FCs = {}
    Ne_ests = {}
    Ne_est_pool = {}
    for Ne_i in Ne_list:
        FCs = [get_FCs(Ne_i, S, n_loci, t, one_subpop, initial_frequencies, no_m) for _ in range(repeats)]
        Ne_FCs[Ne_i] = FCs
        Ne_ests[Ne_i] = [get_ests(FC[0],S, t) for FC in FCs]
        Ne_est_pool[Ne_i] = get_ests(np.mean(np.array(FCs)), S, t)

    ## Varying S
    print("S")
    S_FCs = {}
    S_ests = {}
    S_est_pool = {}
    for S_i in S_list:
        FCs = [get_FCs(Ne, S_i, n_loci, t, one_subpop, initial_frequencies, no_m) for _ in range(repeats)]
        S_FCs[S_i] = FCs
        S_ests[S_i] = [get_ests(FC[0], S_i, t) for FC in FCs]
        S_est_pool[S_i] = get_ests(np.mean(np.array(FCs)), S_i, t)

    # Varying t
    print("t")
    t_FCs = {}
    t_ests = {}
    t_est_pool = {}
    for t_i in t_list:
        FCs = [get_FCs(Ne, S, n_loci, t_i, one_subpop, initial_frequencies, no_m) for _ in range(repeats)]
        t_FCs[t_i] = FCs
        t_ests[t_i] = [get_ests(FC[0], S, t_i) for FC in FCs]
        t_est_pool[t_i] = get_ests(np.mean(np.array(FCs)), S, t_i)


    # varying number of loci
    print("N_loci")
    loc_FCs = {}
    loc_ests = {}
    loc_est_pool = {}
    for loc_i in loc_list:
        FCs = [get_FCs(Ne, S, loc_i, t, one_subpop, initial_frequencies, no_m) for _ in range(repeats)]
        loc_FCs[loc_i] = FCs
        loc_ests[loc_i] = [get_ests(FC[0], S, t) for FC in FCs]
        loc_est_pool[loc_i] = get_ests(np.mean(np.array(FCs)), S, t)


    # Varying initial frequency
    print("Freq")
    freq_FCs = {}
    freq_ests = {}
    freq_est_pool = {}
    for freq_i in freq_list:
        FCs = [get_FCs(Ne, S, n_loci, t, one_subpop, freq_i, no_m) for _ in range(repeats)]
        freq_FCs[freq_i[0]] = FCs
        freq_ests[freq_i[0]] = [get_ests(FC[0], S, t) for FC in FCs]
        freq_est_pool[freq_i[0]] = get_ests(np.mean(np.array(FCs)), S, t)


    # Pooled plotting

    def add_sub_pooled(p_list, est, lab):
        plt.plot(p_list, est, 'k-')
        plt.plot(p_list, [1 for _ in p_list], 'k--')
        plt.xticks(p_list)
        plt.xlabel(lab)
        #plt.ylabel(r'$\hat{N}_e / N_e$')
        plt.ylim(0,1.2)

    fig = plt.figure()
    plt.title("Allelic fluctuation estimator performace")
    plt.subplot(2,5,1)
    Ne_est_pooled_list  = [Ne_est_pool[Ne_i]/Ne_i for Ne_i in Ne_list]
    add_sub_pooled(Ne_list, Ne_est_pooled_list, "Ne")


    plt.subplot(2,5,2)
    S_est_pooled_list = [S_est_pool[S_i]/Ne for S_i in S_list]
    add_sub_pooled(S_list, S_est_pooled_list, "S")


    plt.subplot(2,5,3)
    t_est_pooled_list = [t_est_pool[t_i]/Ne for t_i in t_list]
    add_sub_pooled(t_list, t_est_pooled_list, "t")

    plt.subplot(2,5,4)
    loc_est_pooled_list = [loc_est_pool[loc_i]/Ne for loc_i in loc_list]
    add_sub_pooled(loc_list, loc_est_pooled_list, "number of loci")

    plt.subplot(2,5,5)
    freq_est_pooled_list = [freq_est_pool[freq_i[0]]/Ne for freq_i in freq_list]
    add_sub_pooled([f[0] for f in freq_list], freq_est_pooled_list, "initial frequency")

    # unpooled plotting

    def add_sub_unpooled(p_list, est_array, lab, w):
        plt.violinplot(est_array.T, positions = p_list, widths = w / 1.8)
        plt.plot(p_list, np.mean(est_array, axis =1), 'k-')
        plt.plot(p_list, np.ones_like(p_list), 'k--')
        plt.xticks(p_list)
        plt.xlabel(lab)
        #plt.ylabel(r'$\hat{N}_e/N_e$')
        plt.ylim(max(-1,np.min(est_array)), min(3, np.max(est_array)))


    plt.subplot(2,5,6)
    est_array = np.array([np.array(Ne_ests[Ne_i])/Ne_i for Ne_i in Ne_list])
    add_sub_unpooled(Ne_list, est_array, 'Ne', 30)


    plt.subplot(2,5,7)
    est_array = np.array([S_ests[S] for S in S_list])/Ne
    add_sub_unpooled(S_list, est_array, 'S', 20)


    plt.subplot(2,5,8)
    est_array = np.array([t_ests[t] for t in t_list])/Ne
    add_sub_unpooled(t_list, est_array, 't', 2)


    plt.subplot(2,5,9)
    est_array = np.array([loc_ests[loc_i] for loc_i in loc_list])/Ne
    add_sub_unpooled(loc_list, est_array, 'number of loci', 60)

    plt.subplot(2,5,10)
    est_array = np.array([freq_ests[freq_i[0]] for freq_i in freq_list])/Ne
    add_sub_unpooled([f[0] for f in freq_list], est_array, 'initial frequency', 0.3)

    ax = fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax.set_ylabel(r'$\hat{N}_e / N_e$', labelpad=20)



    plt.show()


    # n_fcs = {}
    # for n in n_subpops:
    #     m_fcs = {}
    #     for m in ms:
    #         print("\nMigration rate = ", m)
    #         m_fcs[m] = []
    #         for r in range(repeats):
    #             print('\t', str(r+1), ' / ', repeats )
    #             FCs = get_FCs(Ne, S, n_loci, gens, n, initial_frequencies, m)
    #             m_fcs[m].append(FCs)
    #         n_fcs[n] = m_fcs


    # maxest = 0
    # for n in n_subpops:
    #     m_fcs = n_fcs[n]
    #     FCmean_Nes = [gens / (2*(np.mean(m_fcs[m]) - 1/(2*S) - 1/(2*S) )) for m in ms] 
    #     plt.plot(ms, FCmean_Nes, '-', label = n)
    #     #plt.scatter(ms, FCmean_Nes, edgecolors='black', color='white')
    #     maxest = max(max(FCmean_Nes), maxest)

    



    # n_fcs = {}
    # for s in S:
    #     m_fcs = {}
    #     for m in ms:
    #         print("\nMigration rate = ", m)
    #         m_fcs[m] = []
    #         for r in range(repeats):
    #             print('\t', str(r+1), ' / ', repeats )
    #             FCs = get_FCs(Ne, s, n_loci, gens, n_subpops, initial_frequencies, m)
    #             m_fcs[m].append(FCs)
    #         n_fcs[s] = m_fcs


    # maxest = 0

    # for s in S:
    #     m_fcs = n_fcs[s]
    #     FCmean_Nes = [gens / (2*(np.mean(m_fcs[m]) - 1/(2*s) - 1/(2*s) )) for m in ms] 
    #     plt.plot(ms, FCmean_Nes, '-', label = s)
    #     #plt.scatter(ms, FCmean_Nes, edgecolors='black', color='white')
    #     maxest = max(max(FCmean_Nes), maxest)

    # plt.plot([min(ms), max(ms)], [Ne,Ne], 'k--')

    # plt.ylim(50,maxest*1.1)
    # plt.xlim(.9*min(ms), 1.1*max(ms))
    # plt.xscale("log")
    # plt.xticks(ticks = ms,labels = ms)
    # plt.xlabel('Migration rate ')
    # plt.ylabel(r'$\hat{N}_e$')
    # plt.legend()

    # plt.show()

    return

if __name__ == "__main__":
    main()
    sys.exit()