import matplotlib.pyplot as plt
import numpy as np

x, simpson, gauss = np.loadtxt("gamma_cdf.txt", unpack=True)

plt.figure(figsize=(10, 6), dpi=100)

plt.plot(x, simpson, '-o', markersize=6, markeredgecolor='black', label='Simpson')
plt.plot(x, gauss, '-<', markersize=6, markeredgecolor='black', label='Gauss 7')

plt.xlabel('x', fontsize=12)
plt.ylabel('Integral Value', fontsize=12)
plt.title("Cumulative Density Function for the Gamma Distribution \n"
          r'$\alpha$'" = 7.5, " r'$\beta$' " = 1.0", 
          pad=20, fontsize=13)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=11, framealpha=1)

plt.tight_layout()
plt.savefig('gammma_results.pdf', bbox_inches='tight')