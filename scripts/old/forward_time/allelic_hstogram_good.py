from allelic import get_FCs
import simuPOP as sim
from simuPOP.sampling import drawRandomSample
from itertools import combinations
import numpy as np





Ne = 100
S = 50
n_loci = 1000
gens = 3
repeats = 100
initial_frequencies = [0.5,0.5]
n_subpops = 1
m = 0
import matplotlib.pyplot as plt


##########################
def fc_variant(freqs_init, freqs_end):
    xi_1, yi_1 = freqs_init[0], freqs_end[0]
    xi_2, yi_2 = freqs_init[1], freqs_end[1]

    return (1/2) * (xi_1-yi_1)**2 / ( (xi_1 + yi_1)/2 - xi_1*yi_1)  + (1/2) * (xi_2-yi_2)**2 / ( (xi_2 + yi_2)/2 - xi_2*yi_2)
def fc_to_Ne(fc, gens, S):
    return gens / (2*(fc - 1/(2*S) - 1/(2*S) ))

##########################

Ne_ests_l = []
FCs = []
for r in range(repeats):
    print(r+1)

    FCs = get_FCs(Ne, S, n_loci, gens, n_subpops, initial_frequencies, m)
    Ne_ests_l.append([fc_to_Ne(FC, gens, S) for FC in FCs])



plt.hist(np.array(Ne_ests_l).flatten(), range = (-Ne,5*Ne), bins = 100)
plt.plot([Ne,Ne], [0,20])
plt.plot([np.mean(np.array(Ne_ests_l)),np.mean(np.array(Ne_ests_l))], [0,20], 'w--')

plt.xlim(-Ne,3*Ne)
plt.show()