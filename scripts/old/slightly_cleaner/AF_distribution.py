from AF_functions_coalescent import get_AF_estimate
import matplotlib.pyplot as plt
import numpy as np


Ne = 100
S = 50
n_subpops=2
m = 0
repeats = 100
t = 3
n_loci = 5000
param_dict = {
        'Ne' : [Ne],
        'S' : [S],
        'n_subpops' : [n_subpops],
        }
e = []
for i in range(repeats):
    e.append(get_AF_estimate(Ne, S, n_subpops, m, t, n_loci))
    print(i)


ep = list(filter(lambda x: x>-100 and x<500, e))

plt.hist(np.array(ep).flatten(), bins =30, color='white',edgecolor='black' )

mean = np.mean(list(filter(lambda x: x>0, e)))

plt.plot([mean,mean], [0,60],'k--', label = 'mean')
plt.plot([Ne,Ne], [0,60],'b--', label = '1')
plt.plot([n_subpops*Ne, Ne*n_subpops], [0,60],'r--', label = 'True')
plt.xlim([-100,500])


plt.legend()
plt.show()