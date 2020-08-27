import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import sys

##########################





def main():
    Ne = 100
    Ne_list = [75,100,150, 200]
    S = 50
    S_list = [20,40,60]
    n_loci = 200
    loc_list = [50, 100, 200, 400]
    t = 3
    t_list = [3,5,7]
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
    Ne_r2s = {}
    Ne_ests = {}
    Ne_est_pool = {}
    for Ne_i in Ne_list:
        print(Ne_i)
        r2s = [get_mean_r2(Ne_i, S, n_loci, t, one_subpop, initial_frequencies, no_m) for _ in range(repeats)]
        Ne_r2s[Ne_i] = r2s
        Ne_ests[Ne_i] = [get_Ne(get_r2_drift(r2[0],S)) for r2 in r2s]
        Ne_est_pool[Ne_i] = get_Ne(get_r2_drift(np.mean(np.array(r2s)),S))


    ## Varying S
    S_r2s = {}
    S_ests = {}
    S_est_pool = {}
    for S_i in S_list:
        print(S_i)
        r2s = [get_mean_r2(Ne, S_i, n_loci, t, one_subpop, initial_frequencies, no_m) for _ in range(repeats)]
        S_r2s[S_i] = r2s
        S_ests[S_i] = [get_Ne(get_r2_drift(r2[0],S_i)) for r2 in r2s]
        S_est_pool[S_i] = get_Ne(get_r2_drift(np.mean(np.array(r2s)),S_i))

    # Varying t
    t_r2s = {}
    t_ests = {}
    t_est_pool = {}
    for t_i in t_list:
        print(t_i)
        r2s = [get_mean_r2(Ne, S, n_loci, t_i, one_subpop, initial_frequencies, no_m) for _ in range(repeats)]
        t_r2s[t_i] = r2s
        t_ests[t_i] = [get_Ne(get_r2_drift(r2[0],S)) for r2 in r2s]
        t_est_pool[t_i] = get_Ne(get_r2_drift(np.mean(np.array(r2s)),S))

    loc_r2s = {}
    loc_ests = {}
    loc_est_pool = {}
    for loc_i in loc_list:
        print(loc_i)
        r2s = [get_mean_r2(Ne, S, loc_i, t, one_subpop, initial_frequencies, no_m) for _ in range(repeats)]
        loc_r2s[loc_i] = r2s
        loc_ests[loc_i] = [get_Ne(get_r2_drift(r2[0],S)) for r2 in r2s]
        loc_est_pool[loc_i] = get_Ne(get_r2_drift(np.mean(np.array(r2s)),S))


    freq_r2s = {}
    freq_ests = {}
    freq_est_pool = {}
    for freq_i in freq_list:
        print(freq_i)
        r2s = [get_mean_r2(Ne, S, n_loci, t, one_subpop, freq_i, no_m) for _ in range(repeats)]
        freq_r2s[freq_i[0]] = r2s
        freq_ests[freq_i[0]] = [get_Ne(get_r2_drift(r2[0],S)) for r2 in r2s]
        freq_est_pool[freq_i[0]] = get_Ne(get_r2_drift(np.mean(np.array(r2s)),S))


    ##############
    # Pooled
    ##############

    def add_sub_pooled(p_list, ests, lab):
        plt.plot(p_list, ests, 'k-')
        plt.plot(p_list, [1 for _ in p_list], 'k--')
        plt.xticks(p_list)
        plt.xlabel(lab)
        #plt.ylabel(r'$\hat{N}_e$')
        plt.ylim(0,1.2)

    fig = plt.figure()
    plt.title("LD estimator performace")

    plt.subplot(2,5,1)
    Ne_est_pooled_list =  [Ne_est_pool[Ne_i]/Ne_i for Ne_i in Ne_list]
    add_sub_pooled(Ne_list, Ne_est_pooled_list, 'Ne')

    plt.subplot(2,5,2)
    S_est_pooled_list = [S_est_pool[S_i]/Ne for S_i in S_list]
    add_sub_pooled(S_list, S_est_pooled_list, 'S')

    plt.subplot(2,5,3)
    t_est_pooled_list = [t_est_pool[t_i]/Ne for t_i in t_list]
    add_sub_pooled(t_list, t_est_pooled_list, 't')

    plt.subplot(2,5,4)
    loc_est_pooled_list = [loc_est_pool[loc_i]/Ne for loc_i in loc_list]
    add_sub_pooled(loc_list, loc_est_pooled_list, 'number of loci')

    plt.subplot(2,5,5)
    freq_est_pooled_list = [freq_est_pool[freq_i[0]]/Ne for freq_i in freq_list]
    add_sub_pooled([f[0] for f in freq_list], freq_est_pooled_list, 'initial frequency')


    ###########
    # UNPOOLED
    ###########

    def add_sub_unpooled(p_list, est_array, lab, d):
        #plt.scatter(np.repeat(S_list, repeats), est_array.ravel(), color='black', s=.5)
        plt.violinplot(est_array.T, positions = p_list, widths = d/1.8)
        plt.plot(p_list, np.mean(est_array, axis =1), 'k-')
        plt.plot(p_list, np.ones_like(p_list), 'k--')
        plt.xlabel(lab)
        #plt.ylabel(r'$\hat{N}_e$')
        plt.ylim(max(0,np.min(est_array)), min(3, np.max(est_array)))
        return


    plt.subplot(2,5,6)
    Ne_ests_array = np.array([np.array(Ne_ests[Ne_i])/Ne_i for Ne_i in Ne_list])
    add_sub_unpooled(Ne_list, Ne_ests_array, 'Ne', 30)

    plt.subplot(2,5,7)
    S_ests_array = np.array([S_ests[S_i] for S_i in S_list])/Ne
    add_sub_unpooled(S_list, S_ests_array, 'S', 20)

    plt.subplot(2,5,8)
    t_ests_array = np.array([t_ests[t_i] for t_i in t_list])/Ne
    add_sub_unpooled(t_list, t_ests_array, 't', 2)

    plt.subplot(2,5,9)
    loc_ests_array = np.array([loc_ests[loc_i] for loc_i in loc_list])/Ne
    add_sub_unpooled(loc_list, loc_ests_array, 'number of loci',60)

    plt.subplot(2,5,10)
    freq_ests_array = np.array([freq_ests[freq_i[0]] for freq_i in freq_list])/Ne
    add_sub_unpooled([f[0] for f in freq_list] , freq_ests_array, 'initial frequency', 0.3)

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