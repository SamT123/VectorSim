#%%
from functions_coalescent import calculate_estimates, get_LD_estimate, get_AF_estimate, stacked_histograms, lines
import matplotlib.pyplot as plt

#%%
Ne = 100
S = 50
n_subpops = 2
n_loci = 1000
t = 3

#%%
params_AF = [Ne, S, n_subpops, t]
params_LD = [Ne, S, n_subpops]
ms = [0.001, 0.01 ,0.1, 0.3,1]
repeats = 50

af_ests = []
ld_ests = [] 

for m in ms:
    print("m = "+str(m))
    print("LD ests")
    ld_ests.append(calculate_estimates(param_list = params_LD, m = m, n_loci = n_loci, repeats = repeats, estimator =  get_LD_estimate))

    print("AF ests")
    af_ests.append(calculate_estimates(param_list = params_AF, m = m, n_loci = n_loci, repeats = repeats, estimator =  get_AF_estimate))



stacked_histograms(af_ests, ms, n_subpops)
plt.show()

stacked_histograms(ld_ests, ms, n_subpops)
plt.show()

lines(af_ests, ld_ests, "AF", "LD", ms)
plt.show()

