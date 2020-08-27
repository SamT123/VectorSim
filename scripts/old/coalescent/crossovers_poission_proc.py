# calculate recombination fractions from recombination rate
# recombination fraction = P(number of crossovers is odd) = P( Pois( lambda ) is odd)
# where lambda = recom_rate * distance between sites
import numpy as np
import matplotlib .pyplot as plt
def get_recom_fraction(rate, distance):
    return 0.5 * (1 - np.exp(- 2 * rate * distance))

distance = 10**np.arange(6,10, .2)
rate = 1e-8


plt.plot(distance,get_recom_fraction(rate, distance))
plt.xscale('log')
plt.show()

recom_fraction_pairs = np.array([0.01,0.1,0.2,0.3,0.4,0.5])
#print(((1-recom_fraction_pairs)**2 + recom_fraction_pairs**2 )/(2*recom_fraction_pairs*( 2 - recom_fraction_pairs)))