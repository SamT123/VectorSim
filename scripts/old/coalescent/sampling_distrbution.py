import island_model
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import KernelDensity

fixed = {'num_replicates' : 10000, 'mutation_rate_scaled' : 1e-5 * 1e6, 'Ne' : 1, 'demes' : 5}

# migration
Ms = np.array([.05])
ms = Ms * (fixed['Ne'] / fixed['demes'])  / (fixed['demes'] - 1)


sims = island_model.run_sim(fixed, ms)
diversities_obs = island_model.get_diversities_all_ms(sims, ms, fixed)

within = diversities_obs[0][0]

true_mean_diversity = np.mean(within)

plt.figure()
for sample_size in [20, 40, 80, 150, 1000]:
    sample_means = np.array([np.mean( np.random.choice(within, size = sample_size, replace = False) ) for i in range(10000)])


    xs = np.arange(0,1010,10)
    ys =  np.exp(KernelDensity(kernel='gaussian', bandwidth = 3).fit(sample_means.reshape(-1, 1)).score_samples(xs.reshape(-1, 1)))
    plt.plot(xs,ys, label = str(sample_size))
#plt.plot([true_mean_diversity]*2, [0,1], color = 'black')
plt.legend()
plt.show()