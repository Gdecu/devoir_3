#ifndef ANALYSE_H
/* ============================================================
 *  Outils d’analyse : convergence – stabilité – énergie – …
 * ============================================================*/
#define ANALYSE_H

#include "model.h"

// Stuff for animation
void anim_enrgy(char *initial_conditions, CSRMatrix *Ksp, CSRMatrix *Msp, double *u, double *v, int n, int T, double dt, int nbr_iter);
void get_coords(double *coord, int n_nodes, char *filename);
void convergence_complexity(CSRMatrix *Ksp, CSRMatrix *Msp, int N, double T, char *init, double *dt_range, int Ndt);



/* ================= stabilité : noyau numérique ============== */
double spectral_radius_newmark(const CSRMatrix *K, const CSRMatrix *M,
                               int n_nodes, double h,
                               double beta, double gamma,
                               int    max_it, double tol);

/* ================= stabilité : ρ(h) (β,γ fixés) ============= */
void stability_scan(const CSRMatrix *K, const CSRMatrix *M,
                    int n_nodes,
                    const double *dt_list, int Ndt,
                    double beta, double gamma);

/* =============== stabilité : ρ(β,γ) (h fixé) ================ */
void stability_map_beta_gamma(const CSRMatrix *K, const CSRMatrix *M,
                              int n_nodes,
                              const double *beta_vec,  int Nbeta,
                              const double *gamma_vec, int Ngamma,
                              double dt_fixed);

/* ==== Petits wrappers de démonstration (facultatifs) : ==== */
void run_stability_vs_dt (const CSRMatrix *K, const CSRMatrix *M,
                          int n_nodes);
void run_stability_vs_bg (const CSRMatrix *K, const CSRMatrix *M,
                          int n_nodes);

#endif   /* ANALYSE_H */