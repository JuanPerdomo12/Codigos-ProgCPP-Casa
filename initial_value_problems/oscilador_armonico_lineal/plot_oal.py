import numpy as np
import matplotlib.pyplot as plt

t1, x1, v1 = np.loadtxt("output-euler.txt", unpack=True)
t2, x2, v2= np.loadtxt("output-heun.txt", unpack=True)

plt.figure(figsize=(10, 8), dpi=100)
plt.plot(x1, v1, '-o', markersize=6, markeredgecolor='black', label='Euler Method')
plt.plot(x2, v2, '-<', markersize=6, markeredgecolor='black', label='Heun Method')

plt.xscale('linear')
plt.yscale('linear')

plt.xlabel('Position (x)', fontsize=12)
plt.ylabel('Velocity (v)', fontsize=12)
plt.title("Phase Diagram of a linear Harmonic Oscillator with " r'$x_0 = 1$' " , " r'$v_0 = 0$' " and " r'$\omega = 4.376$', fontsize=15)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12, framealpha=1)
plt.tight_layout()
plt.savefig("phase-diagram-v-x.pdf", bbox_inches='tight')