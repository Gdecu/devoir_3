#!/usr/bin/env python3
"""
plot_convergence.py
Trace l’erreur ‖u(Δt) − u(ref)‖₂ / ‖u(ref)‖₂
en fonction du pas de temps Δt pour vérifier l’ordre de convergence
de la méthode de Newmark.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Les Δt employés dans la fonction `convergence` :
DT_LIST = [0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0]

# Dossier où `convergence` écrit les fichiers ‘final_dt*.txt’
DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "dt"


def load_displacements(path: Path) -> np.ndarray:
    """
    Lit un fichier 'final_dt*.txt' et renvoie le vecteur des déplacements
    [ux0, uy0, ux1, uy1, …] (on ignore les vitesses).
    On saute les lignes qui ne comportent pas exactement 4 valeurs.
    """
    disp = []
    with path.open() as f:
        for line in f:
            cols = line.split()
            if len(cols) == 4:                      # format : ux uy vx vy
                disp.extend([float(cols[0]), float(cols[1])])
    return np.array(disp)


def main() -> None:
    # Chargement des résultats
    files = [DATA_DIR / f"final_dt{dt:.5g}.txt" for dt in DT_LIST]
    u_all = [load_displacements(fp) for fp in files]

    # solution de référence : Δt = 0.01 (premier fichier)
    u_ref = u_all[0]

    # Calcul de l'erreur relative L2
    errors = [np.linalg.norm(u - u_ref) / np.linalg.norm(u_ref) for u in u_all]

    # Affichage du taux de convergence (pente de la droite log-log)
    coeffs = np.polyfit(np.log(DT_LIST[1:]), np.log(errors[1:]), 1)  # on évite le point de référence
    order = -coeffs[0]
    print(f"Ordre de convergence estimé : {order:.2f}")

    # Trace
    plt.loglog(DT_LIST, errors, "o-", label=f"Pente ≈ {order:.2f}")
    plt.xlabel(r"$\Delta t$")
    plt.ylabel(r"Erreur relative $L^2$")
    plt.title("Convergence temporelle – méthode de Newmark")
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.4)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
