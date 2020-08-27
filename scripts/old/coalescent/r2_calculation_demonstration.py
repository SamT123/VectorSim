 #RECOMBINATION

# length: units are arbitrary and continuous
# recombiantion_rate: X_over / unit / generation
import msprime
import tskit
import numpy as np
tree_sequence = msprime.simulate(
    sample_size = 50, Ne = 100, length = 1e8, recombination_rate = 2e-8, mutation_rate = 1e-10
)

def count_trees(tree_sequence):
    return sum([1 for tree in tree_sequence.trees()])

def count_mutations(tree_sequence):
    return sum([ sum([1 for site in tree.sites()]) for tree in tree_sequence.trees() ])

def myLDcalc(tree_sequence, ):
    for variant in tree_sequence.variants():
        if variant.site.id == 0:
            gen1 = variant.genotypes
        if variant.site.id == 1:
            gen2 = variant.genotypes
            break
    n = len(gen1)
    pA = sum(gen1) / n
    pB = sum(gen2) / n
    pAB= sum(np.logical_and(gen1, gen2))/n
    pab= sum(np.logical_and(1-gen1, 1-gen2))/n
    pAb= sum(np.logical_and(gen1, 1-gen2))/n
    paB= sum(np.logical_and(1-gen1, gen2))/n

    print(sum([pAB, pAb, paB, pab]))

    r_unphased = ( pAB - pA * pB ) / np.sqrt(pA * (1-pA) * pB * (1-pB))
    r2_phased   = (pAB * pab - pAb * paB)**2 / (pA * (1-pA) * pB * (1-pB))
    return r_unphased**2, r2_phased

for variant in tree_sequence.variants():
    print(
        variant.site.id, variant.site.position,
        variant.alleles, sep = '\t'
    )

print(count_mutations(tree_sequence))
print(tskit.LdCalculator(tree_sequence).r2_matrix())
print(tskit.LdCalculator(tree_sequence).r2(0,1))
print(myLDcalc(tree_sequence))