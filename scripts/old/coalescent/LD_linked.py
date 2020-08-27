# E(r2) = 1 / 3Ne + 1 / S wihtout migration
# and we maybe can ignore the sampling term

# actually, get full N estimator from Hill 1981 paper. 
# check for sensititivity of Waples simplification to recombination rate.


# THEN test in multiple deme scenario

# check sensitivity to migration rate
# compare to their analytical expectation for bias this would cause
# is it robust to wider demographic scenrios than they tested?

#######################################################################
import operator
import sys
import msprime
import tskit
import numpy as np
import itertools
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
import matplotlib.pyplot as plt
######################################
# LD to estimate a single population #
######################################


# handling tree_sequences

def print_first(tree_sequence):
    print(tree_sequence.first().draw(format = 'unicode'))
    return

def count_trees(tree_sequence):
    return sum([1 for tree in tree_sequence.trees()])

def count_mutations(tree_sequence):
    return sum([ sum([1 for site in tree.sites()]) for tree in tree_sequence.trees() ])

def get_allele_filter(tree_sequence, alpha,S):
    return np.array([  np.abs(0.5 - sum(v.genotypes)/S) < 0.5 - alpha for v in tree_sequence.variants()])

def get_LD_matrix(tree_sequence):
    return tskit.LdCalculator(tree_sequence).r2_matrix()

def filter_matrix(LD_matrix, freq_filter):
    return LD_matrix[freq_filter][:, freq_filter]

def get_locations(tree_sequence):
    return np.array([v.site.position for v in tree_sequence.variants()])


def get_recom_fraction(rate, distance):
    return 0.5 * (1 - np.exp(- 2 * rate * distance))

def get_pairwise_distances(tree_sequence):
    n_sites = count_mutations(tree_sequence)
    locations = get_locations(tree_sequence)

    pairwise_distances = np.zeros((n_sites, n_sites))

    for i, j in itertools.product(range(n_sites), range(n_sites)):
        d = np.abs(locations[i] - locations[j])
        pairwise_distances[i,j] = d

    return pairwise_distances


def get_pairwise_recom_fractions(tree_sequence, recom_rate ,flat_list = False):
    n_sites = count_mutations(tree_sequence)
    locations = get_locations(tree_sequence)

    fractions = np.zeros((n_sites, n_sites))
    flat = []
    for i, j in itertools.product(range(n_sites), range(n_sites)):
        d = np.abs(locations[i] - locations[j])
        frac = get_recom_fraction(recom_rate, d)
        fractions[i,j] = frac
        if i != j:
            flat.append(frac)

    if flat_list:
        return flat

    return fractions

def remove_outliers(arr, min_val, max_val):
    return arr[np.logical_and(arr > min_val, arr < max_val)]

def summarise_sims(sims):
    sims_mean = np.mean(sims)
    sims_med  = np.median(sims)
    sims_mean_clean = np.mean(remove_outliers(sims, 0,10000))
    sims_sd = np.sqrt(np.var(sims)) / np.sqrt(len(sims))

    print("mean:\t\t" + str(sims_mean) + ' +- ' + str(2*sims_sd) + '\n' + 'median:\t\t' + str(sims_med) + '\n' + "clean mean:\t" + str(sims_mean_clean))




# ESTIMATORS 

# The simple Waples one: Ne_hat = 1 / 3[ E(r2) - 1/S ]
# This assumes all recombination fractions = 0.5

def get_r2_drift(LD_mat,S):
    r2_drift_matrix = (LD_mat - 1/(1*S) ) / (1 - 1/(1*S))
    return r2_drift_matrix

def off_diag_mean(matrix):
    return (np.sum(matrix) - np.trace(matrix)) / ( (matrix.shape[0] - 1) * matrix.shape[1])

def get_off_diag_elems(matrix):
    idx = np.where( np.logical_and(~np.eye(matrix.shape[0],dtype=bool), np.triu(np.ones_like(matrix)) ))
    return matrix[idx]

def Ne_waples(LD_mat):
    E_r2 = off_diag_mean(LD_mat)
    return  1 / ( 3 * (E_r2) )

def get_pairs(tree_sequence):
    return list(itertools.combinations(range(count_mutations(tree_sequence)), 2))

def Ne_hill(tree_sequence):
    global S

    allele_frequency_filter = get_allele_filter(tree_sequence,.1,S)
    accepted = np.where(allele_frequency_filter)[0]


    LD_matrix = get_LD_matrix(tree_sequence)
    recom_fraction_matrix = get_pairwise_recom_fractions(tree_sequence, recom_rate)
    #site_pairs = get_pairs(tree_sequence)
    site_pairs = list(itertools.combinations(accepted,2))


    recom_fraction_pairs = np.zeros(len(site_pairs))
    LD_pairs = np.zeros(len(site_pairs))

    i = 0
    for pair in site_pairs:
        recom_fraction_pairs[i] = recom_fraction_matrix[pair[0], pair[1]]
        LD_pairs[i] = LD_matrix[pair[0], pair[1]]
        i +=1


    #plt.hist(recom_fraction_pairs, bins = 100)
    #plt.show()
    gammas = ((1-recom_fraction_pairs)**2 + recom_fraction_pairs**2 )/(2*recom_fraction_pairs*( 2 - recom_fraction_pairs))
    #Âgammas = np.full_like(gammas,1/3)


    alpha = (LD_pairs - 1/S)/gammas


    var_alphas = 2 * np.mean(LD_pairs)**2 / gammas**2
    #var_alphas = 1

    N_inv = sum(alpha / var_alphas) / ((len(alpha)/var_alphas) ) # sum !!!!!!!!!!!!!!!!!
    N_inv = sum(alpha / var_alphas) / sum(1/var_alphas)  # sum !!!!!!!!!!!!!!!!!



    return 1 / N_inv


#filtered_matrix = filter_matrix(LD_matrix, allele_frequency_filter)





######################################
######################################

# important to consider how to set up coalescent tree
# size of genomic window is important, as it affects recombination fraction
# then we can control mutation rate to decide the number of mutant sites to consider

# chromosome is 100 cM long

def main():

    S                 = 50
    recom_rate        = 1e-8
    mutation_rate     = 2e-9 * 100
    Ne_arr            = np.array([50, 100, 150, 200])
    length            = 1e8
    reps              = 3

    WaplesEsts = np.zeros(shape = (len(Ne_arr), reps) )
    mean_r2s = np.zeros(shape = (len(Ne_arr), reps) )
    mean_r2_drifts = np.zeros(shape = (len(Ne_arr), reps) )
    raw_r2s = []
    raw_physical_dists = []

    for Ne_idx in range(len(Ne_arr)):

        Ne = Ne_arr[Ne_idx]
        print('Running simulations for Ne = ' + str(Ne))

        raw_r2s.append([])
        raw_physical_dists.append([])

        for repeat in range(reps):
            print(repeat + 1)

            tree_sequence = msprime.simulate(sample_size=S, Ne=Ne, length=length, recombination_rate=recom_rate, mutation_rate=mutation_rate/Ne, model = "dtwf") 

            count_trees(tree_sequence)

            allele_frequency_filter = get_allele_filter(tree_sequence, 0.3, S)

            LD_matrix = get_LD_matrix(tree_sequence)
            filtered_matrix = filter_matrix(LD_matrix, allele_frequency_filter)
            filtered_r2_drift = get_r2_drift(filtered_matrix,S)
            raw_r2s[-1].append(get_off_diag_elems(filtered_r2_drift))


            pairwise_distances = get_pairwise_distances(tree_sequence)
            filtered_pairwise_distances = filter_matrix(pairwise_distances, allele_frequency_filter)
            raw_physical_dists[-1].append(get_off_diag_elems(filtered_pairwise_distances))

            WapEst = Ne_waples(filtered_r2_drift)
            mean_r2 = off_diag_mean(filtered_matrix)
            mean_r2_drift = off_diag_mean(filtered_r2_drift)


            WaplesEsts[Ne_idx, repeat] = WapEst
            mean_r2s[Ne_idx, repeat] = mean_r2
            mean_r2_drifts[Ne_idx, repeat] = mean_r2_drift



    np.savetxt("../data/WaplesEsts.csv", WaplesEsts)
    np.savetxt("../data/r2s.csv", mean_r2s )
    np.savetxt("../data/r2_drifts.csv", mean_r2_drifts )

    # print(mean_r2s)
    # plt.scatter( np.repeat(Ne_arr, reps), mean_r2_drifts.flatten(), color = 'black', marker = '_')
    # plt.plot( Ne_arr, np.mean(mean_r2_drifts, axis = 1), 'k-', label = r'$\mathdefault{mean }\:r^2_{drift}$')
    # plt.plot( Ne_arr, np.mean(mean_r2s, axis = 1), 'k--', label = r'$\mathdefault{mean }\,r^2$')
    # plt.plot( Ne_arr, 1 / (3*Ne_arr), 'k:', label = r'$\frac{1}{3N_e}$')
    # plt.xticks(Ne_arr)
    # plt.legend(loc = 'upper right')
    # plt.show()








    plt.scatter( np.repeat(Ne_arr, reps), mean_r2s.flatten(), color = 'black', marker = '_')
    plt.plot( Ne_arr, np.mean(mean_r2s, axis = 1), 'k-', label = r'$\mathrm{mean }\,r^2$')
    plt.plot( Ne_arr, 1/(2*S) + (1 - (1/(2*S)) ) * 1 / (3*Ne_arr), 'k--', label = r'$\frac{1}{3N_e} + E(r^2_{sample})$')
    plt.plot( Ne_arr, 1 / (3*Ne_arr), 'k:', label = r'$\frac{1}{3N_e}$')
    plt.xticks(Ne_arr)
    plt.legend(loc = 'upper right')
    #plt.show()
    lincols = ['aqua', 'yellow', 'pink', 'orange']





    ######################
    # Recombination rate #
    ######################

    for Ne_idx in range(len(raw_r2s)):
        plt.subplot(2,2,Ne_idx+1)
        for repeat_idx in range(reps):


            x = raw_physical_dists[Ne_idx][repeat_idx]
            x2= get_recom_fraction(recom_rate, x)
            y = raw_r2s[Ne_idx][repeat_idx]

            L = sorted(zip(x2,y), key=operator.itemgetter(0))
            new_x, new_y = zip(*L)
            y2 = lowess(new_y, new_x, frac = .03,  return_sorted=False)



            title = 'Ne = ' + str(Ne_arr[Ne_idx])

            plt.scatter(x2, y, marker = '.', s = 4, alpha = .5, color = 'black')
            plt.plot(new_x, y2, color = 'aqua')
            plt.title(title)
            plt.xlabel('physical distance')
            plt.ylabel('r2 value')
            plt.xlim(-0.1,.6)
            #plt.xscale('log')
        c=np.arange(0,.5,.01)
        y_theory_sved = 1/(1+4*Ne_arr[Ne_idx]*c)
        y_theory_hill = ((1-c)**2+ c**2)/(2*Ne_arr[Ne_idx]*c*(2-c)) 
        plt.plot(c, y_theory_hill, color = 'red', label = "Hill")
        plt.plot(c, y_theory_sved, color = 'pink', label = 'Sved')
        plt.legend()


    plt.show()


if __name__ == "__main__":
    main()
    sys.exit()