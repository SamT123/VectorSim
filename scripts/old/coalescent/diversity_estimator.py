import msprime 
import tskit
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import itertools
# sampling strategy and population structure:

# fixed parameters


# sampling from 2 populations only
def get_configs(demes):
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=10),
        msprime.PopulationConfiguration(sample_size=10),
        ] + [ msprime.PopulationConfiguration(sample_size=0) for i in range(0, demes-2)]
    return population_configurations




def run_sim(fixed_params, ms):

    sims = []
    i=0
    for m in ms:

        print('.'*i); i+=1
        population_configurations = get_configs(fixed_params['demes'])

        migration_matrix = np.full((fixed_params['demes'], fixed_params['demes']), m )
        np.fill_diagonal(migration_matrix, 0)
        migration_matrix = migration_matrix.tolist()

        replicates = msprime.simulate(
            Ne = fixed_params['Ne'],
            population_configurations = population_configurations,
            migration_matrix = migration_matrix,
            mutation_rate = fixed_params['mutation_rate_scaled'],
            length = 1,
            num_replicates = fixed_params['num_replicates']
            )
        sims.append(replicates)
    return sims

def get_diversities(sim, fixed_params):
    within_diversity = np.zeros(shape = fixed_params['num_replicates'])
    between_diversity = np.zeros(shape = fixed_params['num_replicates'])

    for i, treeseq in enumerate(sim):
        indicies = [list(range(5*n,5*n+5)) for n in range(0, fixed_params['demes'])]

        
        between_div = treeseq.divergence(sample_sets = [range(0,10), range(10,20)])
        within_div = treeseq.diversity(sample_sets = [range(0,10), range(10,20)])
        

        between_diversity[i] = np.mean(between_div)
        within_diversity[i] = np.mean(within_div)
    return within_diversity, between_diversity

def get_diversities_all_ms(sims, ms, fixed_params):
    within_diversities = np.zeros(shape = (len(ms),fixed_params['num_replicates']))
    between_diversities = np.zeros(shape = (len(ms),fixed_params['num_replicates']))
    i = 0
    for sim in sims:
        
        divs = get_diversities(sim, fixed_params)
        within_diversities[i] = divs[0]
        between_diversities[i] = divs[1]
        print('.'*i); i+=1

    return within_diversities, between_diversities

def get_analyticals(fixed_params, Ms):

    within_expected = 4 * fixed_params['Ne'] * fixed_params['demes']  * fixed_params['mutation_rate_scaled'] * np.ones_like(Ms)
    between_expected = 4 * fixed_params['Ne'] * fixed_params['demes']  * fixed_params['mutation_rate_scaled'] * (1 + 1 / Ms)

    return within_expected, between_expected

def get_error_bars(obs):
    l = []
    m = []
    u = []
    for m_data in obs:
        print(len(m_data))

        mid = np.mean(m_data)
        lower = mid - 1.96 * np.std(m_data) / np.sqrt(len(m_data))
        upper = mid + 1.96 * np.std(m_data) / np.sqrt(len(m_data))
        
        l.append(lower)
        m.append(mid)
        u.append(upper)
    return l,m,u
    
def main():

    fixed = {'num_replicates' : 100, 'mutation_rate_scaled' : 1e-5 * 1e6, 'Ne' : 1, 'demes' : 50}

    # migration
    Ms = np.arange(.05, 0.3, .05)
    ms = Ms   / (fixed['demes'] - 1)



    sims = run_sim(fixed, ms)
    diversities_obs = get_diversities_all_ms(sims, ms, fixed)
    l_w, m_w, u_w = get_error_bars(diversities_obs[0])
    l_b, m_b, u_b = get_error_bars(diversities_obs[1])

    diversities_an  = get_analyticals(fixed, ms)

    print("Complete")

    plt.subplot(1,2,1)
    plt.plot(Ms, m_w, color = 'black', label = 'Observed')
    plt.plot(Ms, l_w, ':', color = 'grey')
    plt.plot(Ms, u_w, ':', color = 'grey')

    plt.plot(Ms, diversities_an[0], '--', color = 'black', label = 'Analytical')
    plt.ylim(0, np.max(np.concatenate([np.mean(diversities_obs[0], axis = 1), diversities_an[0]])*1.1))
    plt.xlabel("mutation rate")
    plt.ylabel("diversity")
    plt.legend()
    plt.title("Within deme genetic diversity")

    plt.subplot(1,2,2)
    plt.plot(Ms, np.mean(diversities_obs[1], axis = 1), color = 'black', label = 'Observed')
    #plt.plot(Ms, diversities_an[1], '--', color = 'black', label = 'Analytical')
    plt.ylim(0, np.max(np.mean(diversities_obs[1], axis = 1)*1.1))

    plt.xlabel("mutation rate")
    plt.legend()

    plt.title("Between deme genetic diversity")

    plt.show()

    return

if __name__ == "__main__":
    main()
    sys.exit()


# print("\nWITHIN POPULATION\n" + "-"*17)
# print("Observed:\t" + str(round(np.mean(within_diversity))))
# print("Analytic:\t" + str(within_expected))


# print("\nbetween POPULATION\n" + "-"*17)
# print("Observed:\t" + str(round(np.mean(between_diversity))))
# print("Analytic:\t" + str(between_expected))

# # histogram of observed diversities
# plt.hist(within_diversity, bins = 30)
# plt.plot(np.full(6, within_diversity_mean), range(0,600,100), 'r--')
# plt.plot(np.full(6, within_expected), range(0,600,100), 'b--')
# #plt.show()

# plt.hist(between_diversity, bins = 30)
# plt.plot(np.full(6, between_diversity_mean), range(0,600,100), 'r--')
# plt.plot(np.full(6, between_expected), range(0,600,100), 'b--')
# #plt.show()


# diversity vs params

# between deme diversities
