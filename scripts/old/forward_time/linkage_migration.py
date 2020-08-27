from linkage_functions import get_migration_matrix, get_mean_r2, get_mean_r2, get_Ne, get_r2_drift, hmean, calculate_estimates
import numpy as np
import matplotlib.pyplot as plt

def make_p_dic(p_list, pooled, ms, fn = lambda x:x):
    p_dic = {k:[] for k in p_list}
    for p in p_list:
        for m in ms:
            p_dic[p].append(fn(pooled[m][p]))
    return p_dic


def add_sub_pooled(p_list, p_dic, lab,ms):
    for p in p_list:
        plt.plot(ms, p_dic[p], label = p)

    plt.plot(ms, np.ones_like(ms), 'k--')
    plt.legend(title = lab)
    plt.xscale('log')
    plt.xticks(ticks = ms, labels= ms)
    #plt.xlabel("Migration rate")
    #plt.ylabel(r'$\hat{N}_e / N_e$')
    return 

def main():
    param_dict = {
        'Ne' : [100],
        'S' : [50],
        'n_loci' : [50],
        't' : [7],
        'n_subpops' : [20],
        'initial_frequencies' : [0.5]
        }

    param_dict_vary = {
    'Ne' : ['Ne',[75,100,150, 200]],
    'S' : ['S',[20,40,60]],
    'n_loci' : ['n_loci', [30,50, 100]],
    't' : ['t', [3,5,7,9]],
    'n_subpops' : ['n_subpops', [10,20,50]],
    'initial_frequencies' : ['initial_frequencies', [0.5,0.75,0.95]]
    }

    ms = [0.001, 0.01 ,0.1,0.3, 49/50]
    repeats = 20
    #print(calculate_estimates(param_dict, ["Ne", [75,100,150, 200]], m))

    m_Ne_ests = {}
    m_Ne_ests_pool = {}
    m_S_ests = {}
    m_S_ests_pool = {}
    m_t_ests = {}
    m_t_ests_pool = {}
    m_n_loci_ests = {}
    m_n_loci_ests_pool = {}
    m_initial_frequencies_ests = {}
    m_initial_frequencies_ests_pool = {}
    m_n_subpops_ests = {}
    m_n_subpops_ests_pool = {}



    for m in ms:
        print("m = ", str(m))
        m_Ne_ests_pool[m], m_Ne_ests[m] = calculate_estimates(param_dict, m,repeats, param_dict_vary['Ne'])
        # m_S_ests_pool[m], m_S_ests[m] = calculate_estimates(param_dict, param_dict_vary['S'], m,repeats)
        # m_t_ests_pool[m], m_t_ests[m] = calculate_estimates(param_dict, param_dict_vary['t'], m,repeats)
        # m_n_loci_ests_pool[m], m_n_loci_ests[m] = calculate_estimates(param_dict, param_dict_vary['n_loci'], m,repeats)
        # m_initial_frequencies_ests_pool[m], m_initial_frequencies_ests[m] = calculate_estimates(param_dict, param_dict_vary['initial_frequencies'], m,repeats)
        # m_n_subpops_ests_pool[m], m_n_subpops_ests[m] = calculate_estimates(param_dict, param_dict_vary['n_subpops'], m,repeats)

    fig = plt.figure()
    plt.subplot(1,2,1)
    Ne_dic_pool = make_p_dic(param_dict_vary['Ne'][1], m_Ne_ests_pool, ms)
    Ne_dic_mean = make_p_dic(param_dict_vary['Ne'][1], m_Ne_ests, ms, fn = np.mean)
    add_sub_pooled(param_dict_vary['Ne'][1], Ne_dic_pool, "Ne",ms)
    plt.subplot(1,2,2)
    add_sub_pooled(param_dict_vary['Ne'][1], Ne_dic_mean, "Ne",ms)
    plt.show()
    return
    plt.subplot(1,6,2)
    S_dic = make_p_dic(S_list, m_S_est_pool)
    add_sub_pooled(S_list, S_dic_corr, "S",ms)

    plt.subplot(1,6,3)
    t_dic = make_p_dic(t_list, m_t_est_pool)
    add_sub_pooled(t_list, t_dic_corr, "t",ms)

    plt.subplot(1,6,4)
    loc_dic = make_p_dic(loc_list, m_loc_est_pool)
    add_sub_pooled(loc_list, loc_dic_corr, r"$N_{loci}$",ms)

    plt.subplot(1,6,5)
    freq_dic = make_p_dic([f[0] for f in freq_list], m_freq_est_pool)
    add_sub_pooled([f[0] for f in freq_list], freq_dic_corr, "Initial Freq",ms)

    plt.subplot(1,6,6)
    ns_dic = make_p_dic(ns_list, m_ns_est_pool)
    add_sub_pooled(ns_list, ns_dic_corr, "Sub-populations",ms)


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
