import matplotlib.pyplot as plt
import numpy as np
# import seaborn as sns
# sns.set_context("poster")

# read data
x, exact, richardson = np.loadtxt("campo.txt", unpack=True)

# plot 
fig, ax = plt.subplots()
ax.plot(x, richardson, '-o', label='Richardson')
ax.plot(x, exact, '-*', label='Exact')
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlabel('Distancia')
ax.set_ylabel('Campo Electrico')
ax.legend()
plt.title(' ')
plt.grid(True)
fig.savefig('campo_vs_derivpotencial.pdf')