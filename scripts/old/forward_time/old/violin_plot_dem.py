import numpy as np
import matplotlib.pyplot as plt

r = np.random.normal(size = (5,50))
x = [1,2,3,5,7]

plt.violinplot(r.T, positions = x)
plt.show()

