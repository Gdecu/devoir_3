#!/usr/bin/env python3
"""
plot_energy.py
Trace l’énergie cinétique, potentielle et totale à partir d’un fichier enrgy.txt
(3 colonnes : Ecin  Epot  Etot)
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent          # dossier src/
FILENAME = HERE.parent / "data" / "enrgy.txt"   # ../data/enrgy.txt
DELIMITER = None        # None => espaces ou tabulations indifféremment
SKIP_HEADER_LINES = 0   # mets 1 si tu as une ligne d’en-tête
# ---------------------------------------------------------------------------


def main(fname: Path):
    if not fname.is_file():
        sys.exit(f"❌ Fichier ‘{fname}’ introuvable")

    # Chargement des 3 colonnes
    try:
        data = np.loadtxt(
            fname,
            delimiter=DELIMITER,
            skiprows=SKIP_HEADER_LINES,
        )
    except ValueError as err:
        sys.exit(f"❌ Problème de format dans {fname} : {err}")

    if data.ndim != 2 or data.shape[1] < 3:
        sys.exit("❌ Le fichier doit comporter 3 colonnes : Ecin  Epot  Etot")

    # Axe des itérations = numéro de ligne (0, 1, 2, ...)
    it   = np.arange(data.shape[0])
    ecin = data[:, 0]
    epot = data[:, 1]
    etot = data[:, 2]

    # Tracé
    plt.figure(figsize=(10, 6))
    plt.plot(it, ecin, label="Énergie Cinétique")
    plt.plot(it, epot, label="Énergie Potentielle")
    plt.plot(it, etot, label="Énergie Totale")

    plt.xlabel("Itération Temporelle")
    plt.ylabel("Énergie")
    plt.title("Évolution de l’Énergie au cours du temps")
    plt.legend()
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main(FILENAME)
