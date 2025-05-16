import numpy as np
import matplotlib.pyplot as plt

# === PARAMÈTRES ===
folder = './data/dt'
"""n_h = 9  # nombre d'étapes de temps
dt_range = [0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1]
T = 20"""
n_h = 17  # nombre d'étapes de temps
dt_range = [0.0001, 0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 5, 10, 15, 20, 50, 100]
T = 200

# === DONNÉES ===
datas = [""]*n_h
for i in range(n_h):
    print(f'Chargement des données pour dt={dt_range[i]}')
    datas[i] = np.loadtxt(f'{folder}/final_dt{dt_range[i]}.txt')

data_ref = datas[0] # On va prendre la solution avec dt=0.01 comme référence

# === ERREURS ===
errors = np.zeros(n_h-1)
for i in range(1, n_h):
    errors[i-1] = np.linalg.norm(datas[i] - data_ref) / np.linalg.norm(data_ref) # Erreur relative
    # errors[i-1] = np.linalg.norm(datas[i] - data_ref) # Erreur absolue
    print(f'Erreur pour dt={dt_range[i]}: {errors[i-1]}')

# === CONVERGENCE ===
# On va tracer l'erreur en fonction de dt
plt.figure()
plt.loglog(dt_range[1:], errors, 'o-')
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('dt')
plt.ylabel('Erreur relative')
plt.title('Convergence de la méthode')
plt.grid()
plt.savefig('convergence.png')
