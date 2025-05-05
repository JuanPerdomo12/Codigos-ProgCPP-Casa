import matplotlib.pyplot as plt
import numpy as np
# import seaborn as sns
# sns.set_context("poster")

# read data
counter, parentesis1, k1, coma, k2, parentesis2, norm  = np.loadtxt("twin_primes.txt", unpack=True)

# plot 
fig, ax = plt.subplots()
ax.plot(counter, norm, '-o', label='')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlabel('nth twin primes')
ax.set_ylabel('norm of twin primes')
ax.legend()
plt.title(' ')
plt.grid(True)
fig.savefig('twin_primes.pdf')