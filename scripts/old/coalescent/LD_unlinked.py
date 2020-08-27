import tskit
import LD_linked as Waples2011
import numpy as np
import matplotlib.pyplot as plt
import msprime 

reps = 10
Ne_arr = [25, 50, 75, 100, 150]
S=100
mutation_rate = 1e-8
recom_rate = 1e-8
positions = [0, 1e8-1, 1e8, 2e8-1]
rates = [recom_rate, 0.5, recom_rate, 0]
num_loci = int(positions[-1])

recombination_map = msprime.RecombinationMap(
    positions=positions, rates=rates, num_loci=num_loci)


r2_observed = []
r2_corrected = []
Ne_pred = []

for Ne_idx in range(len(Ne_arr)):
    Ne = Ne_arr[Ne_idx]
    r2_observed_list = []
    r2_corrected_list = []
    Ne_list = []
    print('Running simulations for Ne = ' + str(Ne))

    for rep in range(reps):
        print(str(rep + 1) + ' / ' + str(reps))
        tree_sequence = msprime.simulate(
            sample_size=S, Ne=Ne, recombination_map=recombination_map,
            model="dtwf", mutation_rate = mutation_rate)

        first_chrom2 = 'not assigned'
        n1=0
        n2=0
        for variant in tree_sequence.variants():
            if variant.position > 1e8:
                first_chrom2 = variant.index
                n2+=1
            if variant.position < (1e8 - 1):
                n1+=1

        r2 = tskit.LdCalculator(tree_sequence).r2_matrix()
        allele_frequency_filter = Waples2011.get_allele_filter(tree_sequence, 0.2, S)
        chrom2_filter = np.logical_and(allele_frequency_filter, np.append(np.zeros(n1),np.ones(n2)))
        chrom1_filter = np.logical_and(allele_frequency_filter, np.append(np.ones(n1),np.zeros(n2)))
        filtered_matrix = r2[chrom1_filter][:, chrom2_filter]

        obs_r2_drift = np.mean(Waples2011.get_r2_drift(filtered_matrix,S))
        obs_r2_uncorrected = np.mean(filtered_matrix)

        N_prediction = 1/(3*obs_r2_drift)

        r2_corrected_list.append(obs_r2_drift)
        r2_observed_list.append(obs_r2_uncorrected)

        Ne_list.append(N_prediction)

    r2_observed.append(np.array(r2_observed_list))
    r2_corrected.append(np.array(r2_corrected_list))

    Ne_pred.append(np.array(Ne_list))


r2_true_Hill = 1 / (3*np.array(Ne_arr))
r2_obs_Hill  = 1/(1*S) + (1 - (1/(1*S)) ) * r2_true_Hill
r2_true_Sved = 1 / (1+2*np.array(Ne_arr))
r2_obs_Sved  = 1/(1*S) + (1 - (1/(1*S)) ) * r2_true_Sved

# r2 values
plt.scatter( np.repeat(Ne_arr, reps), np.array(r2_observed).flatten(), color = 'black', marker = '_')
plt.scatter( np.repeat(Ne_arr, reps), np.array(r2_corrected).flatten(), color = 'blue', marker = '_')

plt.plot( Ne_arr, np.mean(r2_observed, axis = 1), 'k-', label = r'$\mathrm{obs }\,r^2$')
plt.plot( Ne_arr, np.mean(r2_corrected, axis = 1), 'b-', label = r'$\mathrm{adj }\,r^2$')


plt.plot( Ne_arr, r2_true_Hill, 'b--', label = r'$\frac{1}{3N_e}$')
plt.plot( Ne_arr, r2_obs_Hill, 'k--', label = r'$\frac{1}{3N_e} + E(r^2_{sample})$')

#plt.plot( Ne_arr, r2_true_Sved, 'k:', label = r'$\frac{1}{2N_e} + E(r^2_{sample})$')



plt.xticks(Ne_arr)
plt.legend(loc = 'upper right')
plt.show()

# Ne prediction
for i in range(len(Ne_arr)):
    plt.subplot(3,2,i+1)
    plt.hist(Ne_pred[i])
plt.show()


