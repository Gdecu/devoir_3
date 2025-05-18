import numpy as np 
import matplotlib.pyplot as plt

fileloc = './data/results/time.txt'
# Load the data
data = np.loadtxt(fileloc)
# Extract time, ux, uy, vx, vy
time = data[:, 0]
ux = data[:, 1]
uy = data[:, 2]
vx = data[:, 3]
vy = data[:, 4]

time_num = np.diff(time)
vx_num = np.diff(ux) / time_num
vy_num = np.diff(uy) / time_num

fig, axs = plt.subplots(2, 2, figsize=(14, 8), sharex=True)

# En haut à gauche : ux(t) et vx(t)
axs[0, 0].plot(time, ux, label=r'$u_x$', color='blue')
axs[0, 0].plot(time, vx, label=r'$v_x$', color='green')
axs[0, 0].legend()
axs[0, 0].grid()

# En haut à droite : uy(t) et vy(t)
axs[0, 1].plot(time, uy, label=r'$u_y$', color='orange')
axs[0, 1].plot(time, vy, label=r'$v_y$', color='red')
axs[0, 1].legend()
axs[0, 1].grid()

# En bas à gauche : ux(t) et vx_num(t)
axs[1, 0].plot(time, ux, label=r'$u_x$', color='blue')
axs[1, 0].plot(time[1:], vx_num, label=r'$v_{x,num}$', color='purple', linestyle='--')
axs[1, 0].set_xlabel('Time (s)')
axs[1, 0].legend()
axs[1, 0].grid()

# En bas à droite : uy(t) et vy_num(t)
axs[1, 1].plot(time, uy, label=r'$u_y$', color='orange')
axs[1, 1].plot(time[1:], vy_num, label=r'$v_{y,num}$', color='brown', linestyle='--')
axs[1, 1].set_xlabel('Time (s)')
axs[1, 1].legend()
axs[1, 1].grid()

# Calcul des bornes communes
x_min, x_max = np.min(time), np.max(time)
y_min = min(np.min(ux), np.min(uy), np.min(vx), np.min(vy), np.min(vx_num), np.min(vy_num))
y_max = max(np.max(ux), np.max(uy), np.max(vx), np.max(vy), np.max(vx_num), np.max(vy_num))

for ax in axs.flat:
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

plt.tight_layout()
plt.savefig('./images/time_evolution_subplots.png')
#plt.show()
