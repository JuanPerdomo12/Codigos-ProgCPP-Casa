import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt("millikan_data.txt", unpack=True)

with open("data.txt") as f:
    m, delta_m, b, delta_b, h_measured, h_error, phi, phi_error = map(float, f.read().split())

plt.figure(figsize=(10, 6), dpi=100)

plt.plot(x, y, 'bo', markersize=6, markeredgecolor='black', label='Datos experimentales')

x_fit = np.linspace(min(x), max(x), 500)
y_fit = m * x_fit + b
plt.plot(x_fit, y_fit, 'r-', linewidth=2, label=f'Ajuste lineal: V = {m:.3e}ν {b:.3f}')

plt.xlabel('Frecuencia ν (Hz)', fontsize=12)
plt.ylabel('Voltaje V (voltios)', fontsize=12)
plt.title("Experimento del Efecto Fotoeléctrico de Millikan \n"
          f"h = ({h_measured:.3e} ± {h_error:.1e}) J·s   "
          f"ϕ = ({phi:.2f} ± {phi_error:.2f}) V", 
          pad=20, fontsize=13)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=11, framealpha=1)

plt.tight_layout()
plt.savefig('millikan_results.pdf', bbox_inches='tight')
plt.show()