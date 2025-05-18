#!/usr/bin/env python3
# ==============================================================
#   STABILITY PLOTS for Newmark scheme
#       • ρ(h)  for β = ¼, γ = ½      (figure 1)
#       • ρ(β,γ) heat-map             (figure 2)
# ==============================================================

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.grid(True, which='both', ls='--', alpha=0.5)
# ------------------------------------------------------------------
#  FIGURE 1 : ρ vs time-step h  (β = ¼, γ = ½)
# ------------------------------------------------------------------
try:
    h, rho = np.loadtxt("./data/spectral_radius_dt.txt", unpack=True)
except OSError:
    print("❌  spectral_radius_dt.txt manquant – lancez run_stability_vs_dt() d’abord")
    h = rho = np.array([])

if h.size:
    fig1, ax1 = plt.subplots()
    ax1.plot(h, rho, "o-", label=r"$\rho(A(h))$")
    ax1.axhline(1.0, ls="--", lw=1, label=r"$\rho=1$ (stable)")
    ax1.set_xscale("log")
    ax1.minorticks_on()                        # graduations mineures
    ax1.grid(True, which="both", ls="--", alpha=0.3)  # grille majeure + mineure
    ax1.set_axisbelow(True)                   # pour que la grille passe sous les courbes

    ax1.set_xlabel(r"pas de temps $h$")
    ax1.set_ylabel(r"rayon spectral $\rho$")
    ax1.set_title(r"Stabilité pour $(\beta,\gamma)=(\frac{1}{4},\frac{1}{2})$")
    ax1.legend()
    ax1.grid(alpha=.3)
    fig1.tight_layout()
    fig1.savefig("./images/rho_vs_dt.png", dpi=300)
    print("✅  rho_vs_dt.png généré")

# ------------------------------------------------------------------
#  FIGURE 2 : heat-map ρ(β,γ)   (h fixé)
# ------------------------------------------------------------------
try:
    beta, gamma, rho_bg = np.loadtxt("./data/spectral_radius_bg.txt", unpack=True)
except OSError:
    print("❌  spectral_radius_bg.txt manquant – lancez run_stability_vs_bg() d’abord")
    beta = gamma = rho_bg = np.array([])

if beta.size:
    B_uni = np.unique(beta)
    G_uni = np.unique(gamma)
    NB, NG = len(B_uni), len(G_uni)
    assert NB * NG == beta.size, "Le fichier n’est pas sur une grille régulière"

    Beta   = beta.reshape(NB, NG)
    Gamma  = gamma.reshape(NB, NG)
    Rho    = rho_bg.reshape(NB, NG)

    fig2, ax2 = plt.subplots(figsize=(6,4))
    pcm = ax2.pcolormesh(Beta, Gamma, Rho, shading="auto")
    # ------------------------------------------------------------
    #  Droite  β = 14 (γ + 1/2)^2  (pour γ ≥ 1/2)
    # ------------------------------------------------------------
    gamma_line = np.linspace(0.5, Gamma.max(), 200)       # γ ≥ 1/2
    beta_line  = 1/4 * (gamma_line + 0.5)**2

    # On ne garde que la partie qui tombe dans la fenêtre du graphique
    inside = (beta_line >= Beta.min()) & (beta_line <= Beta.max())
    ax2.plot(beta_line[inside], gamma_line[inside],
            "w--", lw=1.5, label=r"$\beta = \frac{1}{4}\,(\gamma+\frac{1}{2})^2$")
    cb  = fig2.colorbar(pcm, ax=ax2, label=r"$\rho$")

    #condition théorique de stabilité
    
    # couple symplectique (¼,½)
    ax2.plot(0.25, 0.5, "r*", ms=10,
             label=r"$(\frac{1}{4},\frac{1}{2})$")

    ax2.set_xlabel(r"$\beta$")
    ax2.set_ylabel(r"$\gamma$")
    ax2.set_title(r"Carte du rayon spectral $\rho$ (h = 0.1)")
    ax2.legend()
    ax2.grid(alpha=.3)
    fig2.tight_layout()
    fig2.savefig("./images/rho_beta_gamma.png", dpi=300)
    print("✅  rho_beta_gamma.png généré")


plt.show()
