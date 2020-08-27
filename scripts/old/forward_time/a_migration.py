from allelic import *
import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import sys


def main():

    Ne = 100
    Ne_list = [75,100,150, 200]
    S = 50
    S_list = [20,40,60]
    n_loci = 50
    loc_list = [30,50, 100]
    t = 3
    t_list = [3,5,7,9]
    repeats = 20
    n_subpops = 20
    ns_list = [10,20,50]
    one_subpop = 1
    ms = [0.001, 0.01 ,0.1,0.3]
    no_m = 0
    initial_frequencies = [0.5,0.5]
    freq_list = [[0.5,0.5],[0.75,0.25], [0.95,0.05]]


    m_Ne_ests = m_Ne_est_pool = {}
    m_S_ests = m_S_est_pool = {}
    m_t_ests = m_t_est_pool = {}
    m_loc_ests = m_loc_est_pool = {}
    m_freq_ests = m_freq_est_pool = {}
    m_ns_ests = m_ns_est_pool = {}

    for m in ms:
        print("m = ", str(m))

        print("\tNe simutlations")
        Ne_FCs = {}
        Ne_ests = {}
        Ne_est_pool = {}
        for Ne_i in Ne_list:
            FCs = [get_FCs(Ne_i, S, n_loci, t, n_subpops, initial_frequencies, m) for _ in range(repeats)]
            Ne_FCs[Ne_i] = FCs
            Ne_ests[Ne_i] = [get_ests(FC[0],S, t) for FC in FCs]
            Ne_est_pool[Ne_i] = get_ests(np.mean(np.array(FCs)), S, t)
        m_Ne_ests[m] = Ne_ests
        m_Ne_est_pool[m] = Ne_est_pool


        ## Varying S
        print("\tS simutlations")
        S_FCs = {}
        S_ests = {}
        S_est_pool = {}
        for S_i in S_list:
            FCs = [get_FCs(Ne, S_i, n_loci, t, n_subpops, initial_frequencies, m) for _ in range(repeats)]
            S_FCs[S_i] = FCs
            S_ests[S_i] = [get_ests(FC[0], S_i, t) for FC in FCs]
            S_est_pool[S_i] = get_ests(np.mean(np.array(FCs)), S_i, t)
        m_S_ests[m] = S_ests
        m_S_est_pool[m] = S_est_pool

        # Varying t
        print("\tt simutlations")
        t_FCs = {}
        t_ests = {}
        t_est_pool = {}
        for t_i in t_list:
            FCs = [get_FCs(Ne, S, n_loci, t_i, n_subpops, initial_frequencies, m) for _ in range(repeats)]
            t_FCs[t_i] = FCs
            t_ests[t_i] = [get_ests(FC[0], S, t_i) for FC in FCs]
            t_est_pool[t_i] = get_ests(np.mean(np.array(FCs)), S, t_i)

        m_t_ests[m] = t_ests
        m_t_est_pool[m] = t_est_pool

        print("\tloci simutlations")
        loc_FCs = {}
        loc_ests = {}
        loc_est_pool = {}
        for loc_i in loc_list:
            FCs = [get_FCs(Ne, S, loc_i, t, n_subpops, initial_frequencies, m) for _ in range(repeats)]
            loc_FCs[loc_i] = FCs
            loc_ests[loc_i] = [get_ests(FC[0], S, t) for FC in FCs]
            loc_est_pool[loc_i] = get_ests(np.mean(np.array(FCs)), S, t)
        m_loc_ests[m] = loc_ests
        m_loc_est_pool[m] = loc_est_pool

        print("\tFrequency simutlations")
        freq_FCs = {}
        freq_ests = {}
        freq_est_pool = {}
        for freq_i in freq_list:
            FCs = [get_FCs(Ne, S, n_loci, t, n_subpops, freq_i, m) for _ in range(repeats)]
            freq_FCs[freq_i[0]] = FCs
            freq_ests[freq_i[0]] = [get_ests(FC[0], S, t) for FC in FCs]
            freq_est_pool[freq_i[0]] = get_ests(np.mean(np.array(FCs)), S, t)
        m_freq_ests[m] =freq_ests
        m_freq_est_pool[m] = freq_est_pool

        print("\tSubpopulations simutlations")
        ns_FCs = {}
        ns_ests = {}
        ns_est_pool = {}
        for ns_i in ns_list:
            FCs = [get_FCs(Ne, S, n_loci, t, ns_i, initial_frequencies, m) for _ in range(repeats)]
            ns_FCs[ns_i] = FCs
            ns_ests[ns_i] = [get_ests(FC[0], S, t) for FC in FCs]
            ns_est_pool[ns_i] = get_ests(np.mean(np.array(FCs)), S, t)
        m_ns_ests[m] = ns_ests
        m_ns_est_pool[m] = ns_est_pool

    def make_p_dic(p_list, pooled):
        p_dic = {k:[] for k in p_list}
        for p in p_list:
            for m in ms:
                p_dic[p].append(pooled[m][p])
        return p_dic

    def add_sub_pooled(p_list, p_dic, lab):
        for p in p_list:
            plt.plot(ms, p_dic[p], label = p)

        plt.plot(ms, np.ones_like(ms), 'k--')
        plt.legend(title = lab)
        plt.xscale('log')
        plt.xticks(ticks = ms, labels= ms)
        #plt.xlabel("Migration rate")
        #plt.ylabel(r'$\hat{N}_e / N_e$')
        
        return 

    fig = plt.figure()
    plt.title("Allelic Fluctuations")
    plt.subplot(1,6,1)
    Ne_dic = make_p_dic(Ne_list, m_Ne_est_pool)
    Ne_dic_corr = {k:np.array(v)/k for k,v in Ne_dic.items()}
    add_sub_pooled(Ne_list, Ne_dic_corr, "Ne")

    plt.subplot(1,6,2)
    S_dic = make_p_dic(S_list, m_S_est_pool)
    S_dic_corr = {k:np.array(v)/Ne for k,v in S_dic.items()}
    add_sub_pooled(S_list, S_dic_corr, "S")

    plt.subplot(1,6,3)
    t_dic = make_p_dic(t_list, m_t_est_pool)
    t_dic_corr = {k:np.array(v)/Ne for k,v in t_dic.items()}
    add_sub_pooled(t_list, t_dic_corr, "t")

    plt.subplot(1,6,4)
    loc_dic = make_p_dic(loc_list, m_loc_est_pool)
    loc_dic_corr = {k:np.array(v)/Ne for k,v in loc_dic.items()}
    add_sub_pooled(loc_list, loc_dic_corr, r"$N_{loci}$")

    plt.subplot(1,6,5)
    freq_dic = make_p_dic([f[0] for f in freq_list], m_freq_est_pool)
    freq_dic_corr = {k:np.array(v)/Ne for k,v in freq_dic.items()}
    add_sub_pooled([f[0] for f in freq_list], freq_dic_corr, "Initial Freq")

    plt.subplot(1,6,6)
    ns_dic = make_p_dic(ns_list, m_ns_est_pool)
    ns_dic_corr = {k:np.array(v)/Ne for k,v in ns_dic.items()}
    add_sub_pooled(ns_list, ns_dic_corr, "Sub-populations")


    ax = fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax.set_ylabel(r'$\hat{N}_e / N_e$', labelpad=20)

    ax = fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax.set_xlabel('Migration rate', labelpad=20)


    plt.show()









    return


if __name__ == "__main__":
    main()
