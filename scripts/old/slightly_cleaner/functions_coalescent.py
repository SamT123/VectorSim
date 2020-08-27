import itertools
import numpy as np
import msprime
import tskit
import itertools
import matplotlib.pyplot as plt

def fc_variant(genotype, S_hap):
    # allele 1 :
    xi_1 = sum(genotype[:S_hap])/S_hap
    yi_1 = sum(genotype[S_hap:])/S_hap

    xi_2 = 1-xi_1
    yi_2 = 1-yi_1

    if xi_1 > .3 and xi_1 < .7:
        return (1/2) * (xi_1-yi_1)**2 / ( (xi_1 + yi_1)/2 - xi_1*yi_1)  + (1/2) * (xi_2-yi_2)**2 / ( (xi_2 + yi_2)/2 - xi_2*yi_2)
    return "nope"

def simulate_snp(Ne, S, num_snps, n_subpops,m, t):

    population_configurations = [msprime.PopulationConfiguration(sample_size=None, initial_size = Ne) for i in range(n_subpops)]
    samples = [ msprime.Sample(population=0, time=0) for i in range(S)] + [msprime.Sample(population=0, time=t) for i in range(S)]
    M = get_migration_matrix(m, n_subpops)

    replicates = msprime.simulate(
        samples=samples, Ne=Ne, mutation_rate=1,
        num_replicates=100 * num_snps, population_configurations=population_configurations, migration_matrix=M)
    t_max = 0
    variants = np.empty((num_snps, 2*S), dtype="u1")
    total_branch_length = np.empty(num_snps)
    j = 0
    num_adaptive_updates = 0
    num_rejected_trees = 0
    for ts in replicates:
        tree = next(ts.trees())
        tbl = tree.get_total_branch_length()
        if tbl > t_max:
            new_t_max = tbl
            new_variants = np.empty((num_snps, 2*S), dtype="u1")
            new_total_branch_length = np.empty(num_snps)
            keep = np.where(np.random.random(j) < t_max / new_t_max)[0]
            j = keep.shape[0]
            new_variants[:j] = variants[keep]
            new_total_branch_length[:j] = total_branch_length[keep]
            variants = new_variants
            total_branch_length = new_total_branch_length
            t_max = new_t_max
            num_adaptive_updates += 1
        else:
            if np.random.random() < tbl / t_max:
                total_branch_length[j] = tbl
                for variant in ts.variants():
                    variants[j] = variant.genotypes
                    break
                else:
                    continue
                    raise Exception("Must have at least one mutation")
                j += 1
                if j == num_snps:
                    break
            else:
                num_rejected_trees += 1
    assert j == num_snps
    #print("num adaptive updates: ", num_adaptive_updates)
    #print("num rejected trees", num_rejected_trees)
    return variants

def get_AF_estimate(Ne, S, n_subpops, t, m, n_loci, rate = None):

    snps = simulate_snp(Ne, S, n_loci, n_subpops, m, t)
    #print(("estimating size...\n"))
    fc_total = 0
    pass_count = 0
    fail_count=0
    for v in snps:
        fc_v = fc_variant(v, S)
        if fc_v != 'nope':
            fc_total += fc_v
            pass_count += 1
        else: fail_count += 1

    #print(pass_count + fail_count)
    fc = fc_total / pass_count
    Ne_est = (t) / (2*(fc - 1/(1*S) - 1/(1*S) ))

    return Ne_est,rate



def calculate_estimates(param_list, m, n_loci, repeats, estimator):
    """ 'param_list' is assumed to be in the same order as the arguments for 'estimator' """
    r2s = []
    ests = []
    rate = False
    for i in range(repeats):
        print(str(i+1) + " / " + str(repeats))
        this_Ne, rate = estimator(*param_list, m = m, rate = rate, n_loci = n_loci)
        ests.append(this_Ne / param_list[0])
    return ests


def get_migration_matrix(m,n):

    if n == 1:
        return np.array([0])
    m_adj = m / (n)
    M = np.full( (n,n), m_adj )
    np.fill_diagonal(M, 0)
    M = M.tolist()
    return M


def get_LD_estimate(Ne, S, n_subpops, m, rate = False, n_loci = 100):

    # set up population parameters

    ## migration matrix
    M = get_migration_matrix(m, n_subpops)
    ## sample
    population_configurations = [msprime.PopulationConfiguration(sample_size=S)] + [msprime.PopulationConfiguration(sample_size=0) for i in range(n_subpops - 1)]
    ## mutation and recombination : adjust based on Ne to get same number of mutations!
    mutation_rate = 0*5e-9 / 40
    recom_rate = 1e-8
    positions = [0, 1e8-1, 1e8, 2e8-1]
    rates = [recom_rate, 0.5, recom_rate, 0]
    num_loci = int(positions[-1])

    recombination_map = msprime.RecombinationMap(
        positions=positions, rates=rates, num_loci=num_loci)

    tree_sequence = msprime.simulate(
        Ne=Ne, recombination_map=recombination_map,
        mutation_rate = mutation_rate, 
        population_configurations=population_configurations,
        migration_matrix=M,
        model = "dtwf")

    if not rate:
        print("Calculating rate...")
        L = 0
        for tree in tree_sequence.trees():
            L += tree.get_length() * tree.get_total_branch_length()
        rate = n_loci/L

    tree_sequence = msprime.mutate(tree_sequence, rate=rate, random_seed=None, model=None, keep=False, start_time=None, end_time=None)

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
    fil = np.zeros((n1+n2,n1+n2))
    fil[:n1, n1:] = 1

    i=0
    for v in tree_sequence.variants():
        if sum(v.genotypes / len(v.genotypes)) < 0.05:
            fil[:,i] = 0
            fil[i,:] = 0
        i += 1


    fil = fil.astype(int)
    r2 = np.where(fil, r2, np.zeros_like(r2))

    r2_mean = np.sum(r2) / np.sum(fil)
    r2_drift = r2_mean -  (1/(1*S) / (1 - 1/(1*S)))
    Ne_est = 1/(3*r2_drift)

    #print(n1+n2)
    #print(Ne_est)
    return Ne_est, rate

def central_mean(d):
    d = list(d)
    d2 = sorted(d)[round(.1*len(d))+1:round(.9*len(d))-1]
    return np.mean(d2)

def lines(dat1, dat2,lab1, lab2, ms):
    dat1_means, dat2_means = list(map(np.mean, dat1)), list(map(np.mean, dat2))
    dat1_centrals, dat2_centrals = list(map(central_mean, dat1)), list(map(central_mean, dat2))

    plt.subplot(2,1,1)
    plt.plot(ms, dat1_means, 'r-', label = lab1)
    plt.plot(ms, dat2_means, 'b-', label = lab2)
    plt.title("Means")
    plt.xscale("log")

    plt.subplot(2,1,2)
    plt.plot(ms, dat1_centrals, 'r-', label = lab1)
    plt.plot(ms, dat2_centrals, 'b-', label = lab2)
    plt.title("Central means")
    plt.xscale("log")
    plt.legend()

    return 


def stacked_histograms(full_dat,ms,n_s):
    rows = len(full_dat)
    i=0
    plt.figure(figsize = [10,15])
    for i in range(len(ms)):
        d = full_dat[i]
        plt.subplot(rows, 1, i+1)
        plt.hist(full_dat[i], bins = np.arange(-.5,5.5,.1), range = [-.5,5.5],color='white',edgecolor='black' )
        mean = central_mean(full_dat[i])
        plt.plot([1,1], [0,50], 'b--')
        plt.plot([1*n_s,1*n_s], [0,50], 'g--')
        plt.plot([mean,mean], [0,50], 'k--')
        plt.title('m = ' + str(ms[i]))
        i+=1

    return