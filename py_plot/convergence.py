import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress

# === PARAMÈTRES ===
folder = './data/convergence'


n_h = 23  # nombre d'étapes de temps
dt_range = [0.0001, 0.001, 0.002, 0.0025, 0.005, 0.0075, 0.01 , 0.015, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 5, 10, 15, 20, 50, 100]
T = 200
n = 192

# === DONNÉES ===
datas = [""]*n_h
for i in range(n_h):
    print(f'Chargement des données pour dt={dt_range[i]}')
    datas[i] = np.loadtxt(f'{folder}/final_dt{dt_range[i]}.txt')
datas = np.array(datas)

print(datas.shape)


data_ref = datas[0] # On va prendre la solution avec dt=0.0001 comme référence
ux_ref = data_ref[:,0]
uy_ref = data_ref[:,1]
vx_ref = data_ref[:,2]
vy_ref = data_ref[:,3]

datas = datas[1:] # On va ignorer la première solution (dt=0.0001) pour la convergence
dt_range = dt_range[1:] # On va ignorer la première solution (dt=0.0001) pour la convergence
n_h -= 1
ux = datas[:,:,0]
uy = datas[:,:,1]
vx = datas[:,:,2]
vy = datas[:,:,3]

# === ERREURS ===
errors_ux = np.zeros((n_h, n))
errors_uy = np.zeros((n_h, n))
errors_vx = np.zeros((n_h, n))
errors_vy = np.zeros((n_h, n))

for i in range(n_h):
    # Erreur relative
    errors_ux[i,:] = abs(ux[i,:] - ux_ref)
    errors_uy[i,:] = abs(uy[i,:] - uy_ref)
    errors_vx[i,:] = abs(vx[i,:] - vx_ref)
    errors_vy[i,:] = abs(vy[i,:] - vy_ref)

errors_ux /= abs(ux_ref)
errors_uy /= abs(uy_ref)
errors_vx /= abs(vx_ref)
errors_vy /= abs(vy_ref)

def regress(dt, err):
    dt = np.log10(dt)
    err = np.log10(err)
    slope, intercept, r_value, p_value, std_err = linregress(dt, err)
    print(f"Ordre estimé de convergence : {slope:.3f}")
    print(f"Constante de proportionnalité : {intercept:.3f}")
    return (10 ** intercept) * (dt ** slope)

# === CONVERGENCE ===
# On va tracer l'erreur en fonction de dt
dt_range = np.array(dt_range) / T

# Choix des noeuds
X = 7
lst = np.random.choice(np.arange(0, n), size=X, replace=False)
print(f"Indices choisis : {lst}")


# Régression linéaire
# Erreur moyenne sur les n = 192 points
mean_error_ux = np.mean(errors_ux, axis=1)
p, C, r_value, p_value, std_err = linregress(np.log10(dt_range[:16]), np.log10(mean_error_ux[:16]))
print("Regression linéaire ux")
print(f"Ordre estimé de convergence : {p:.3f}")
print(f"Constante de proportionnalité : {C:.3f}")
regr_ux =  (10**C) *  (( dt_range[:16] )**(p))

mean_error_uy = np.mean(errors_uy, axis=1)
p, C, r_value, p_value, std_err = linregress(np.log10(dt_range[:16]), np.log10(mean_error_uy[:16]))
print("Regression linéaire uy")
print(f"Ordre estimé de convergence : {p:.3f}")
print(f"Constante de proportionnalité : {C:.3f}")
regr_uy =  (10**C) *  (( dt_range[:16] )**(p))

fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

# --- ux ---
x = 0
for i in lst:
    if x == 0:
        axs[0].loglog(dt_range, errors_ux[:,i], color='gray', alpha=0.4, linewidth=1, marker='o', markersize=2, linestyle='-', label = "Erreur relative avec les noeuds {lst}".format(lst=lst))
        x = 1
    else:
        axs[0].loglog(dt_range, errors_ux[:,i], color='gray', alpha=0.4, linewidth=1, marker='o', markersize=2, linestyle='-')
axs[0].loglog(dt_range, mean_error_ux, label='Erreur moyenne', color='royalblue', linewidth=2.5, marker='o', markersize=5)
axs[0].loglog(dt_range[:16], regr_ux, label='régression', linestyle='--', color='black')
axs[0].set_title(r"Convergence $u_x$")
axs[0].set_xlabel(r"$dt/T = N^{-1}$ (log)")
axs[0].set_ylabel('Erreur relative (log)')
axs[0].grid(True, which='both', linestyle='--', alpha=0.6)
axs[0].legend()

# --- uy ---
x = 0
for i in lst:
    if x == 0:
        axs[1].loglog(dt_range, errors_uy[:,i], color='gray', alpha=0.4, linewidth=1, marker='o', markersize=2, linestyle='-', label = "Erreur relative avec les noeuds {lst}".format(lst=lst))
        x = 1
    else:
        axs[1].loglog(dt_range, errors_uy[:,i], color='gray', alpha=0.4, linewidth=1, marker='o', markersize=2, linestyle='-')
axs[1].loglog(dt_range, mean_error_uy, label='Erreur moyenne', color='darkorange', linewidth=2.5, marker='o', markersize=5)
axs[1].loglog(dt_range[:16], regr_uy, label='régression', linestyle='--', color='black')
axs[1].set_title(r"Convergence $u_y$")
axs[1].set_xlabel(r"$dt/T = N^{-1}$ (log)")
axs[1].grid(True, which='both', linestyle='--', alpha=0.6)
axs[1].legend()

plt.tight_layout()
plt.savefig('./images/convergence_ux_uy.png')
#plt.show()


# --- vx et vy ---

# Erreur moyenne sur les 7 points sélectionnés
mean_error_vx = np.mean(errors_vx[:,lst], axis=1)
mean_error_vy = np.mean(errors_vy[:,lst], axis=1)

# Régression linéaire vx
p, C, r_value, p_value, std_err = linregress(np.log10(dt_range[:16]), np.log10(mean_error_vx[:16]))
print("Regression linéaire vx")
print(f"Ordre estimé de convergence : {p:.3f}")
print(f"Constante de proportionnalité : {C:.3f}")
regr_vx = (10**C) * (dt_range[:16]**p)

# Régression linéaire vy
p, C, r_value, p_value, std_err = linregress(np.log10(dt_range[:16]), np.log10(mean_error_vy[:16]))
print("Regression linéaire vy")
print(f"Ordre estimé de convergence : {p:.3f}")
print(f"Constante de proportionnalité : {C:.3f}")
regr_vy = (10**C) * (dt_range[:16]**p)

fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

# --- vx ---
for i in lst:
    axs[0].loglog(dt_range, errors_vx[:,i], color='gray', alpha=0.4, linewidth=1, marker='o', markersize=2, linestyle='-')
axs[0].loglog(dt_range, mean_error_vx, label='Erreur moyenne', color='seagreen', linewidth=2.5, marker='o', markersize=5)
axs[0].loglog(dt_range[:16], regr_vx, label='régression', linestyle='--', color='black')
axs[0].set_title(r"Convergence $v_x$")
axs[0].set_xlabel(r"$dt/T = N^{-1}$ (log)")
axs[0].set_ylabel('Erreur relative (log)')
axs[0].grid(True, which='both', linestyle='--', alpha=0.6)
axs[0].legend()

# --- vy ---
for i in lst:
    axs[1].loglog(dt_range, errors_vy[:,i], color='gray', alpha=0.4, linewidth=1, marker='o', markersize=2, linestyle='-')
axs[1].loglog(dt_range, mean_error_vy, label='Erreur moyenne', color='firebrick', linewidth=2.5, marker='o', markersize=5)
axs[1].loglog(dt_range[:16], regr_vy, label='régression', linestyle='--', color='black')
axs[1].set_title(r"Convergence $v_y$")
axs[1].set_xlabel(r"$dt/T = N^{-1}$ (log)")
axs[1].grid(True, which='both', linestyle='--', alpha=0.6)
axs[1].legend()

plt.tight_layout()
plt.savefig('./images/onvergence_vx_vy.png')
#plt.show()