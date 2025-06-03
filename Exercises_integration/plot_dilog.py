import matplotlib.pyplot as plt
import numpy as np

x, dilog = np.loadtxt("dilogarithm.txt", unpack=True)

plt.figure(figsize=(10, 6), dpi=100)

plt.plot(x, dilog, '-o', markersize=6, markeredgecolor='black', label='Dilogarithm value')

plt.xlabel('x', fontsize=12)
plt.ylabel('Di-logarithm Value', fontsize=12)
plt.title("Di-Logarithm Function \n", 
          pad=20, fontsize=13)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=11, framealpha=1)

plt.tight_layout()
plt.savefig('dilog_results.pdf', bbox_inches='tight')  