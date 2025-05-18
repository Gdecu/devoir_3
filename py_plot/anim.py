import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# === PARAMÈTRES ===
folder = './data/anim'
n_frames = 201  # nombre d'étapes de temps
coord_file = './data/coords.txt'

# === Chargement des positions initiales ===
initial_coords = np.loadtxt(coord_file, delimiter=",")[:, :2]
n_nodes = initial_coords.shape[0]

# === Chargement des déplacements ===
ux_uy_all = []

for i in range(n_frames):
    filepath = os.path.join(folder, f'final_{i}.txt')
    data = np.loadtxt(filepath)
    ux_uy_all.append(data[:, :2])  # colonnes ux, uy
    print(f"Frame {i+1}/{n_frames} loaded.")

ux_uy_all = np.array(ux_uy_all)  # shape: (n_frames, n_nodes, 2)

# === Positions totales : positions initiales + déplacements ===
amplifFact = 100  # facteur d'amplification des déplacements
positions = ux_uy_all * amplifFact + initial_coords[None, :, :]  # broadcasting

# === Préparation de l'animation ===
fig, ax = plt.subplots()
sc = ax.scatter(positions[0, :, 0], positions[0, :, 1], c='blue', s=20)
ax.set_title("Animation des déplacements absolus")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Fixe les limites de l'axe selon les positions totales
ax.set_xlim(-0.1, 0.1)
ax.set_ylim(-0.08, 0.12)

def update(frame):
    sc.set_offsets(positions[frame])
    return sc,

ani = animation.FuncAnimation(fig, update, frames=n_frames, interval=50, blit=True)

# === Sauvegarde en GIF ===
ani.save("./images/animation.gif", writer='pillow')
