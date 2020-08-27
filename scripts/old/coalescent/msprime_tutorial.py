import msprime
import numpy as np
# BASICS
# make and print tree

tree_sequence = msprime.simulate(sample_size = 6, Ne = 1000)
tree = tree_sequence.first()
print(tree.draw(format = 'unicode'))

# display tree times = travese tree
# each node is an integer, access node_info with tree methods
# root node has parent 'msprime.tskit.NULL'

## eg1. get node times
u = 1
while u != msprime.tskit.NULL:
    print('node {}: time = {}'.format(u, tree.time(u)))
    u = tree.parent(u)

## eg2. get branch length to parent
print(tree.branch_length(2))
print(tree.time(2))
print(tree.time(tree.parent(2)))
print(tree.total_branch_length)

# RECOMBINATION

# length: units are arbitrary and continuous
# recombiantion_rate: X_over / unit / generation

tree_sequence = msprime.simulate(
    sample_size = 6, Ne = 1000, length = 1e4, recombination_rate = 2e-8
)

# multiple trees for different regions of the genome
# iterate over trees with .trees() method - do not store tree objects separately
for tree in tree_sequence.trees():
    print('-'*20)
    print('tree {}: interval = {}'.format(tree.index, tree.interval))
    print(tree.draw(format = 'unicode'))



# VARIANTS
# specify mutation rate and random seed
tree_sequence = msprime.simulate(
    sample_size = 20, Ne = 1e4, length = 5e3, recombination_rate = 2e-8,
    mutation_rate = 2e-8, random_seed = 10
)

# do not need to store /work with entire genomes : just get variant info
for variant in tree_sequence.variants():
    print(
        variant.site.id, variant.site.position,
        variant.alleles, variant.genotypes, sep = '\t'
    )

# HISTORICAL SAMPLES
# can get samples from time in past

## use samples object:
samples = [
    msprime.Sample(population=0, time=0),
    msprime.Sample(0,0),
    msprime.Sample(0,1),
    msprime.Sample(0,1)
] # do I need to construct a big list if there are eg 10,000 samples?

# leave default param values
treeseq = msprime.simulate(samples=samples)
tree = treeseq.first()
print(tree.draw(format = 'unicode'))



# REPLICATION
# checking analytical results --> independent replicates of same simulation
# num_replicates, then simple iteration
# eg lets check the number of segregating sites
num_replicates = 100
n = 10
theta = 5

S = np.zeros(num_replicates)
reps = msprime.simulate(
    Ne = .5,
    sample_size =  10,
    mutation_rate = theta / 2,
    num_replicates = num_replicates
    )


for j, tree_sequence in enumerate(reps):
    S[j] = tree_sequence.num_sites

S_mean_a = np.sum(1/np.arange(1,n )) * theta
S_var_a = (
    theta * np.sum(1 / np.arange(1, n)) +
    theta**2 * np.sum(1 / np.arange(1, n)**2))

print("              mean              variance")
print("Observed      {}\t\t{}".format(np.mean(S), np.var(S)))
print("Analytical    {:.5f}\t\t{:.5f}".format(S_mean_a, S_var_a))


# POPULATION STRUCTURE

# M is the overall symmetric migration rate, and d is the number
# of subpopulations.
num_replicates = 1000000000000
M = 0.2
d = 3
m = M / (2 * (d - 1))
# Allocate the initial sample. Because we are interested in the
# between-subpopulation coalescence times, we choose one sample each
# from the first two subpopulations.
population_configurations = [
    msprime.PopulationConfiguration(sample_size=1),
    msprime.PopulationConfiguration(sample_size=1),
    msprime.PopulationConfiguration(sample_size=0)]
# Now we set up the migration matrix. Since this is a symmetric
# island model, we have the same rate of migration between all
# pairs of subpopulations. Diagonal elements must be zero.
migration_matrix = [
    [0, m, m],
    [m, 0, m],
    [m, m, 0]]
# We pass these values to the simulate function, and ask it
# to run the required number of replicates.
replicates = msprime.simulate(Ne=0.5,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,
    num_replicates=num_replicates)
# And then iterate over these replicates
T = np.zeros(num_replicates)
for i, tree_sequence in enumerate(replicates):
    tree = tree_sequence.first()
    print(tree.time(tree.root))
    T[i] = tree.time(tree.root) / 4
    break
# Finally, calculate the analytical expectation and print
# out the results
analytical = d / 4 + (d - 1) / (4 * M)
print("Observed  =", np.mean(T))
print("Predicted =", analytical)


