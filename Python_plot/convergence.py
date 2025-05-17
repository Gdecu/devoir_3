import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress

# === PARAMÈTRES ===
folder = './data/convergence'
"""n_h = 9  # nombre d'étapes de temps
dt_range = [0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1]
T = 20"""
n_h = 23  # nombre d'étapes de temps
dt_range = [0.0001, 0.001, 0.002, 0.0025, 0.005, 0.0075, 0.01 , 0.015, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 5, 10, 15, 20, 50, 100]
T = 200

# === DONNÉES ===
datas = [""]*n_h
for i in range(n_h):
    print(f'Chargement des données pour dt={dt_range[i]}')
    datas[i] = np.loadtxt(f'{folder}/final_dt{dt_range[i]}.txt')

data_ref = datas[0] # On va prendre la solution avec dt=0.01 comme référence
#print(data_ref[:,:2])
#print(data_ref[:,2:])

# === ERREURS ===
errors_u = np.zeros(n_h-1)
errors_v = np.zeros(n_h-1)
errors = np.zeros(n_h-1)
for i in range(1, n_h):
    errors_u[i-1] = np.linalg.norm(datas[i][:,:2] - data_ref[:,:2]) / np.linalg.norm(data_ref[:,:2]) # Erreur relative sur le déplacement
    errors_v[i-1] = np.linalg.norm(datas[i][:,2:] - data_ref[:,2:]) / np.linalg.norm(data_ref[:,2:]) # Erreur relative sur la vitesse
    errors[i-1] = np.linalg.norm(datas[i] - data_ref) # Erreur absolue
    #print(f'Erreur pour dt={dt_range[i]}: {errors[i-1]}')


# === CONVERGENCE ===
sns.set_theme(style="darkgrid")

# On va tracer l'erreur en fonction de dt
dt_range = np.array(dt_range[1:]) / T

# Supposons que dt_vals et err_vals contiennent les données numériques
# Extraire la zone où dt est entre 1e-4 et 1e-2 :
linear_dt = dt_range[:18]
linear_errors_u = np.log10((errors_u[:18]))
linear_errors_v = np.log10(errors_v[:18])

# Régression linéaire
slope, intercept, r_value, p_value, std_err = linregress(linear_dt, linear_errors_u)
p = np.log10(slope)
p = 1.45
print(f"Ordre estimé de convergence : {p:.3f}")
C = 4
regr = 10**C *( dt_range )**(p)
plt.figure()
plt.loglog(dt_range, errors_u, 'o-', label='Erreur u', markersize=4)
plt.loglog(dt_range, errors_v, 'o-', label='Erreur v', markersize=4)
plt.loglog(dt_range, errors, 'o-', label='Erreur totale', markersize=4)
plt.loglog(dt_range[:16], regr[:16], linestyle='--')
plt.legend()
plt.xlabel(r"Rapport $dt/T = N^{-1}$ (log)")
plt.ylabel('Erreur relative (log)')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.savefig('convergence.png')
