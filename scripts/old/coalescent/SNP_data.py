
import msprime
import numpy as np 
def simulate_snp(sample_size, num_snps):
    replicates = msprime.simulate(
        sample_size=sample_size, Ne=0.5, mutation_rate=10,
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

test=simulate_snp(5,1000)
print(test)
print(np.histogram(np.sum(test,axis=1),bins=range(1,6),normed=True))
print()


theory = 1/np.arange(1,6.)[:-1]
theory /= sum(theory)
print(theory)
