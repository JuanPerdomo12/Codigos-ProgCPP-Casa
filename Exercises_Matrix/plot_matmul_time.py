# plot_matmul_time.py
import numpy as np
import matplotlib.pyplot as plt

N, times = np.loadtxt("matmul_time.txt", unpack=True)

# Normalize times by the first one T(4)
times /= times[0]

plt.figure(figsize=(10, 7), dpi=100)
plt.plot(N, times, '-o', markersize=6, markeredgecolor='black', label='Matrix Multiplication')
theory = (N/N[0])**3 
plt.plot(N, theory, '--', label='Theoretical O(N³)')

plt.xscale('log')
plt.yscale('log')

plt.xlabel('Matrix Size N', fontsize=12)
plt.ylabel(r'Normalized Time $T(N)/T(4)$', fontsize=12)
plt.title("Elapsed Time of Matrix Multiplication", fontsize=13)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=11, framealpha=1)
plt.tight_layout()
plt.savefig("matmul_times.pdf", bbox_inches='tight')
