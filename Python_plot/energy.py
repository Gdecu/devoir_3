import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("./data/enrgy.txt")
dt = 0.1  # Time step
T = 100  # Total time
time = np.arange(0, T, dt)
E_kin = data[:, 0]
E_pot = data[:, 1]
E_tot = data[:, 2]

plt.plot(time, E_kin, label="Kinetic Energy", color="deepskyblue")
plt.plot(time, E_pot, label="Potential Energy", color="orange")
plt.plot(time, E_tot, label="Total Energy", color="green", linestyle="--")
plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("Energy [J]")
plt.title("Conservation of Energy")
plt.grid(True)
plt.savefig("Energy.png")
