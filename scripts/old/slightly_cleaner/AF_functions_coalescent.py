from linkage_functions_coalescent import get_migration_matrix
import msprime


def fc_variant(genotype, S_hap):
    # allele 1 :
    xi_1 = sum(genotype[:S_hap])/S_hap
    yi_1 = sum(genotype[S_hap:])/S_hap

    xi_2 = 1-xi_1
    yi_2 = 1-yi_1

    if xi_1 > .3 and xi_1 < .7:
        return (1/2) * (xi_1-yi_1)**2 / ( (xi_1 + yi_1)/2 - xi_1*yi_1)  + (1/2) * (xi_2-yi_2)**2 / ( (xi_2 + yi_2)/2 - xi_2*yi_2)
    return "nope"



def get_AF_estimate(Ne, S, n_subpops, m, t, n_loci):



    population_configurations = [msprime.PopulationConfiguration(sample_size=None, initial_size = Ne) for i in range(n_subpops)]
    samples = [ msprime.Sample(population=0, time=0) for i in range(S)] + [msprime.Sample(population=0, time=t) for i in range(S)]
    M = get_migration_matrix(m, n_subpops)

    tree_seq = msprime.simulate(length = 2e8, samples=samples, population_configurations = population_configurations, Ne = Ne, mutation_rate = 0, recombination_rate = 2e-9, migration_matrix = M, model = 'dtwf')

    #print(("adding mutations..."))
    L = 0
    for tree in tree_seq.trees():
        L += tree.get_length() * tree.get_total_branch_length()

    tree_seq = msprime.mutate(tree_seq, rate=n_loci/L, random_seed=None, model=None, keep=False, start_time=None, end_time=None)

    #print(("estimating size...\n"))
    fc_total = 0
    count = 0
    bc=0
    for v in tree_seq.variants():
        fc_v = fc_variant(v.genotypes,S)
        if fc_v != 'nope':
            fc_total += fc_v
            count+=1
        else: bc+=1


    fc = fc_total / count
    Ne_est = (t) / (2*(fc - 1/(1*S) - 1/(1*S) ))

    return Ne_est


get_AF_estimate(200,100, 2, 0, 3, 100)

# calculate Fc

# ests = []
# Fcs = []
# Fcs_first = []

# for rep in range(reps):

#     print(rep + 1)
#     tree_seq = msprime.simulate(length = 2e8, samples=samples, Ne = Ne, mutation_rate = mutation_rate, recombination_rate = recom_rate)

#     Fc_total = 0
#     count = 0

#     for v in tree_seq.variants():
#         fcv = fc_variant(v.genotypes,S_hap)
#         if fcv != 'nope':
#             Fc_total += fcv
#             count+=1
#             fcf = fcv

#     Fc = Fc_total / count
#     Fcs.append(Fc)
#     Fcs_first.append(fcf)
#     Ne_est = (t) / (2*(Fc - 1/(2*S_dip) - 1/(2*S_dip) ))

#     ests.append(Ne_est)