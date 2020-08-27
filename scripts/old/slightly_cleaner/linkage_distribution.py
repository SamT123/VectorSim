import matplotlib.pyplot as plt
import numpy as np
from linkage_functions_coalescent import calculate_estimates

Ne = 100
S = 50
n_subpops=2
m = 1
repeats = 100

param_dict = {
        'Ne' : [Ne],
        'S' : [S],
        'n_subpops' : [n_subpops],
        }


e = np.array(calculate_estimates(param_dict, m = m , repeats = repeats, vary = False))
e = e[e<10]
e = e[e>0]
plt.hist(np.array(e).flatten(), bins =30, color='white',edgecolor='black' )
plt.plot([np.mean(e),np.mean(e)], [0,60],'k--', label = 'mean')
plt.plot([1,1], [0,60],'b--', label = '1')
plt.plot([n_subpops, n_subpops], [0,60],'r--', label = 'True')
#plt.xlim([0,10])


plt.legend()
plt.show()