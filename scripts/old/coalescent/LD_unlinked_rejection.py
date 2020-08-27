import msprime
import numpy as np
import itertools
import matplotlib.pyplot as plt
def simulate_snp(Ne, sample_size, num_snps):
    replicates = msprime.simulate(
        sample_size=sample_size, Ne=Ne, mutation_rate=1,
        num_replicates=100 * num_snps)
    t_max = 0
    variants = np.empty((num_snps, sample_size), dtype="u1")
    total_branch_length = np.empty(num_snps)
    j = 0
    num_adaptive_updates = 0
    num_rejected_trees = 0
    for ts in replicates:
        tree = next(ts.trees())
        tbl = tree.get_total_branch_length()
        if tbl > t_max:
            new_t_max = tbl
            new_variants = np.empty((num_snps, sample_size), dtype="u1")
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
    print("num adaptive updates: ", num_adaptive_updates)
    print("num rejected trees", num_rejected_trees)
    return variants


def get_mean_r2(snps):
    S = len(snps[0])
    n_loci = len(snps)
    frequencies = np.sum(snps, axis = 1) / S
    pairs = itertools.combinations(range(n_loci - 1), r = 2)

    r2_tot, count = 0, 0

    for pair in pairs:
        
        pAB = sum(np.logical_and(snps[pair[0]], snps[pair[1]])) / S
        pA, pB = frequencies[pair[0]], frequencies[pair[1]]
        D_i = pAB - pA * pB
        r2_i = D_i**2 / ( pA * (1-pA) * pB * (1-pB) )
        r2_tot += r2_i
        count += 1
    return r2_tot / count


S_hap=50
S_dip= int(S_hap/2)
n_loci = 100
Ne = 10000
repeats = 10
estimates = []
r2s = []

for r in range(repeats):
    print(r+1)
    print("Simulating")
    snps = simulate_snp(Ne, 2*S_dip, n_loci)
    print("r2 calc")
    r2_total = get_mean_r2(snps)
    r2_drift = ( r2_total - 1/(2*S_dip) ) / (1 - 1/(2*S_dip))
    Ne_est   = 1/(3*r2_drift)
    estimates.append(Ne_est)
    r2s.append(r2_total)

print(np.mean(estimates))
print(estimates)
plt.hist(estimates, range = [-Ne, 3*Ne], bins = 300)
plt.plot([Ne,Ne], [0,20])
plt.xlim(-Ne,3*Ne)
plt.show()

print(r2s)
r2_exp = 1/(2*S_dip) + (1-1/(2*S_dip))*(1/(3*Ne))
plt.hist(r2s, bins = 10)
plt.plot([r2_exp,r2_exp], [0,20])
plt.xlim(0, .1)
plt.show()
