from linkage import get_migration_matrix, get_mean_r2, get_mean_r2, get_Ne, get_r2_drift
import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import sys


Ne = 100
Ne_list = [75,100,150, 200]
S = 50
S_list = [20,40,60]
n_loci = 50
loc_list = [30,50, 100]
t = 7
t_list = [3,5,7,9]
repeats = 20
n_subpops = 20
ns_list = [10,20,50]
one_subpop = 1
ms = [0.001, 0.01 ,0.1,0.3]
no_m = 0
initial_frequencies = [0.5,0.5]
freq_list = [[0.5,0.5],[0.75,0.25], [0.95,0.05]]


ests_over_m = {}
exp_over_m = {}
for m in ms:
    r2s = [get_mean_r2(Ne, S, n_loci, t, n_subpops, initial_frequencies, m) for _ in range(repeats)]
    ests_over_m[m] = [get_Ne(get_r2_drift(r2[0],S)) for r2 in r2s]
    exp_over_m[m] = get_Ne( ((1 - m)**2 - m**2/(n_subpops-1) ) / (3*Ne) )

#print(np.mean(ests_over_m[m]))
plt.plot(ms, [np.mean(ests_over_m[m]) for m in ms], 'k-')
plt.plot(ms, [exp_over_m[m] for m in ms], 'k--')
plt.show()


    