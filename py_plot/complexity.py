import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


# === COMPLEXITÉ ===
# On va tracer le temps d'exécution en fonction de dt

# Charger les données de complexité
fileloc = './data/complexity.txt'
complexity_data = np.loadtxt(fileloc)
T = 200
dt_vals = complexity_data[:, 0] 
times = complexity_data[:, 1]

N = T / dt_vals 

# Régression linéaire
slope, intercept, r_value, p_value, std_err = linregress(np.log10(N), np.log10(times))
p = slope
print(f"Ordre estimé de convergence : {p:.3f}")
C = intercept
print(f"Constante de proportionnalité : {C:.3f}")
regr = (10**C) *( N )**(p)

plt.figure()
plt.loglog(N, times, 'o-', label='Temps de calcul', markersize=4)
plt.loglog(N, regr, label='Droite approché', linestyle='--')
plt.xlabel('Nombre d\'iterations N = T / dt (log)')
plt.ylabel('Temps de calcul log(s)')
plt.title('Complexité temporelle en fonction de dt')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.savefig('./images/complexity.png')
