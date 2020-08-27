"""See distribution of estiates for Waples LD estimator under specific single set of parameter values."""

from linkage_functions import get_migration_matrix, get_mean_r2, get_mean_r2, get_Ne, get_r2_drift, hmean, calculate_estimates
import numpy as np
import matplotlib.pyplot as plt

def main():
    param_dict = {
        'Ne' : [100],
        'S' : [50],
        'n_loci' : [100],
        't' : [20],
        'n_subpops' : [2],
        'initial_frequencies' : [0.5]
        }


    m = 1 # panmixia
    repeats = 50

    e_pool, e = calculate_estimates(param_dict, m = m , repeats = repeats, vary = False)
    plt.hist(np.array(e).flatten(), bins =30, color='white',edgecolor='black' )
    plt.plot([np.mean(e),np.mean(e)], [0,60],'k--', label = 'mean')
    #plt.plot([hmean(e),hmean(e)], [0,60],'b--', label = 'harmonic')
    plt.plot([e_pool,e_pool], [0,60],'g--', label = 'pooled')
    plt.legend()
    plt.show()

    return 

if __name__ == "__main__":
    main()