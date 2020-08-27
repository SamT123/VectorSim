import msprime
import numpy as np
import matplotlib.pyplot as plt

# number of generations
t=3
# number of chromosomes sampled
S_dip = 50
S_hap = S_dip * 2
Ne = 200
mutation_rate = 1e-9
length = 2e8
recom_rate = 0
reps = 300

samples = [ msprime.Sample(population=0, time=0) for i in range(S_hap)] + [msprime.Sample(population=0, time=3) for i in range(S_hap)]


def fc_variant(genotype, S_hap):
    # allele 1 :
    xi_1 = sum(v.genotypes[:S_hap])/S_hap
    yi_1 = sum(v.genotypes[S_hap:])/S_hap

    xi_2 = 1-xi_1
    yi_2 = 1-yi_1

    if xi_1 > 40/400 and xi_1 < 360/400:
        if xi_2 > -40/400:
            #print(xi_1,yi_1,xi_2, yi_2)

            return (1/2) * (xi_1-yi_1)**2 / ( (xi_1 + yi_1)/2 - xi_1*yi_1)  + (1/2) * (xi_2-yi_2)**2 / ( (xi_2 + yi_2)/2 - xi_2*yi_2)

    return "nope"


# calculate Fc

ests = []
Fcs = []
Fcs_first = []

for rep in range(reps):

    print(rep + 1)
    tree_seq = msprime.simulate(length = 2e8, samples=samples, Ne = Ne, mutation_rate = mutation_rate, recombination_rate = recom_rate)

    Fc_total = 0
    count = 0

    for v in tree_seq.variants():
        fcv = fc_variant(v.genotypes,S_hap)
        if fcv != 'nope':
            Fc_total += fcv
            count+=1
            fcf = fcv

    Fc = Fc_total / count
    Fcs.append(Fc)
    Fcs_first.append(fcf)
    Ne_est = (t) / (2*(Fc - 1/(2*S_dip) - 1/(2*S_dip) ))

    ests.append(Ne_est)


mean_Fc = sum(Fcs) / len(Fcs)
mean_Fcf = sum(Fcs_first) / len(Fcs_first)

print("mean Fc value:")
print(mean_Fc)
print("\n\nNe_estimate based on mean Fc")
print((t) / (2*(mean_Fc - 1/(2*S_dip) - 1/(2*S_dip) )))
print("\n\nNe_estimate based on first Fc")

print((t) / (2*(mean_Fcf - 1/(2*S_dip) - 1/(2*S_dip) )))




plt.hist(ests,range=(-100,1000), bins = 300, color = 'black')
plt.xlabel(r'$\hat{N}_e$')
plt.ylabel("Frequency")
plt.plot([100,100], [0,70], 'r--')
plt.show()