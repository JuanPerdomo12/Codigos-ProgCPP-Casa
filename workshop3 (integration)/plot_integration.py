import matplotlib.pyplot as plt
import numpy as np
# import seaborn as sns
# sns.set_context("poster")

# read data
N, tr, et, rtr, ert, sr, es, rsr, ers = np.loadtxt("integration.txt", unpack=True)

# plot 
fig, ax = plt.subplots()
ax.plot(N, et, '-o', label='Error trapezoid')
ax.plot(N, ert, '-s', label='Error richardson trapezoid')
ax.plot(N, es, '-<', label='Error simpson')
ax.plot(N, ers, '-*', label='Error richardson simpson')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('$N$')
ax.set_ylabel('Error')
ax.legend()
plt.grid(True)
plt.title('Integration Methods Error in function of N')
fig.savefig('integration.pdf')