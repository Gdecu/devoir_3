#include "devoir_2.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "model.h"
#include "analyse.h"
#include "utils.h"
#include "utils_gmsh.h"
#include "gmshc.h"
#include <time.h>

void newmark_analyse(
    double *u, double *v,          // u et v initiales
    char *intit_cond,
    CSRMatrix *K, CSRMatrix *M,    // CSR K et vecteur M
    double nbr_iter, double dt, int n_nodes
) {
    int N    = 2 * n_nodes;                // DOF
    double beta = 0.25, gamma = 0.5;


    double *p    = malloc(N * sizeof(double));
    double *rhs  = malloc(N * sizeof(double));
    double *q_moy = malloc(N * sizeof(double));
    double *q_new = malloc(N * sizeof(double));
    double *Kq   = malloc(N * sizeof(double));
    double *Aeff = malloc(K->nnz * sizeof(double));            // car M strictement inclu dans K --> nnz K = nnz A
    double *Beff = malloc(K->nnz * sizeof(double));

    //test_build_M_coeffK();        // c bon il marche

    double coeff = beta * dt * dt;                          // coeff pour K coeff : beta h²
    build_M_coeffK(N, K, M, Aeff, coeff);                   // Aeff = ( M + beta h² K )
    
    coeff = - 0.5 * dt * dt * (1.0 - 2.0*beta);             // coeff pour K coeff : - 0.5 h² (1-2beta)
    build_M_coeffK(N, K, M, Beff, coeff);                   // Beff = (M - 0.5 dt² (1-2beta) K )

    double *Ax = malloc(N * sizeof(double));


    FILE *energy_file = fopen("./data/enrgy.txt", "w");
    if (!energy_file) {
        perror("Failed to open energy.txt");
    }

    double kinetic, potential, total;

    for (int k = 0; k < nbr_iter; k++) {
        //printf("iteration %d / %d\n", k+1, nbr_iter);

        //
        // One iteration of the Newmark method
        //

        // 1) p_n = M·v_n
        Matvec(N, M->row_ptr, M->col_idx, M->data, v, p);

        // 2) construction du second membre
        //    rhs = (M - ½h²(1-2beta)·K) q + h·p,

        Matvec(N, K->row_ptr, K->col_idx, Beff, u, rhs);    // rhs = Beff·q = (M - 0.5h²(1-2beta)·K )·q
        cblas_daxpy(N, dt, p, 1, rhs, 1);                   // rhs = rhs + h·p

        // 3) résolution Aeff · q_{n+1} = rhs
        CG(N, K->row_ptr, K->col_idx, Aeff, rhs, q_new, 1e-15);

        // 4) p_{n+1} = p_n - h·K·[ (1-gamma)q_n + gamma q_{n+1} ]
        for (int i = 0; i < N; i++)
            q_moy[i] = (1.0 - gamma) * u[i] + gamma * q_new[i];

        Matvec(N, K->row_ptr, K->col_idx, K->data, q_moy, Kq);
        
        for (int i = 0; i < N; i++)
            p[i] -= dt * Kq[i];

        // 5) mise à jour finale : q ← q_{n+1}, M v ← p 
        cblas_dcopy(N, q_new, 1, u, 1); // q = qnp1
        CG(N, M->row_ptr, M->col_idx, M->data, p, v, 1e-15);


        // 
        // End of one iteration of the Newmark method
        //

        // Calcul de l'énergie cinétique et potentielle - et stockage dans <energy.txt>
        kinetic = 0.0; potential = 0.0; total = 0.0;

        // Ekin = 1/2 v^T M v
        //printf("v[0], v[1] = %le, %le\n", v[0], v[1]);
        //printf("M[0], M[1] = %le, %le\n", M->data[0], M->data[1]);
        Matvec(N, M->row_ptr, M->col_idx, M->data, v, Ax);
        //printf("Ax[0 et 1] = %le, %le\n", Ax[0], Ax[1]);
        kinetic = 0.5 * cblas_ddot(N, v, 1, Ax, 1);

        // Epot = 1/2 u^T K u
        //printf("u[0], u[1] = %le, %le\n", u[0], u[1]);
        //printf("K[0], K[1] = %le, %le\n", K->data[0], K->data[1]);
        Matvec(N, K->row_ptr, K->col_idx, K->data, u, Ax);
        //printf("Ax[0 et 1] = %le, %le\n", Ax[0], Ax[1]);
        potential = 0.5 * cblas_ddot(N, u, 1, Ax, 1);

        // Etot = Ekin + Epot
        total = kinetic + potential;
        fprintf(energy_file, "%.15le %.15le %.15le\n", kinetic, potential, total);
        //printf("Ekin = %.15le, Epot = %.15le, Etot = %.15le\n", kinetic, potential, total);

        // Stockage solution dans <final.txt> : le déplacement et la vitesse au temps T, au même format que le fichier <initial.txt>
        char filename[40];
        sprintf(filename, "./data/anim/final_%d.txt", k+1);
        stock_final(n_nodes, filename, u, v);
    }
    printf("E_total = %.15le\n", total);
    fclose(energy_file);
    free(Ax);

    free(p); free(rhs); free(q_new); free(Kq); free(q_moy); free(Beff);
    free(Aeff);
}

void anim_enrgy(char *initial_conditions, CSRMatrix *Ksp, CSRMatrix *Msp, double *u, double *v, int n, int T, double dt, int nbr_iter) {

    // On stocke l'état intiale
    get_intial_condition(initial_conditions, u, v, 2*n);
    //printf("ux[0]  = %15le, uy[0] = %15le\n", u[0], u[1]);
    //printf("vx[0]  = %15le, vy[0] = %15le\n", v[0], v[1]);
    stock_final(n, "./data/anim/final_0.txt", u, v);


    printf("n  = %d\n", n);
    printf("T  = %d\n", T);
    printf("dt  = %le\n", dt);
    printf("nbr_iter  = %d\n\n", nbr_iter);

    newmark_analyse(u, v,
        initial_conditions,
        Ksp, Msp, 
        nbr_iter, dt , n);
}


void get_coords(double *coord, int n_nodes, char *filename) {

    FILE *coord_file = fopen(filename, "w");
    if (!coord_file) {
        perror("Failed to open coords.txt");
    }

    for (int i = 0; i < n_nodes; i++) {
        fprintf(coord_file, "%.15le, %.15le\n", coord[2 * i], coord[2 * i + 1]);
    }

    fclose(coord_file);
}



void convergence_complexity(CSRMatrix *Ksp, CSRMatrix *Msp, int N, double T, char *init, double *dt_range, int Ndt) {

    double dt = dt_range[0];
    int n_m;
    double *u = (double *)malloc( N * sizeof(double));
    double *v = (double *)malloc( N * sizeof(double));
    double *u_moy = (double *)malloc( N * sizeof(double));
    double *v_moy = (double *)malloc( N * sizeof(double));
    for (int i = 0; i < N; i++) {
        u_moy[i] = 0.0;
        v_moy[i] = 0.0;
    }
    double *t = (double *)malloc( (((int) T / dt ) + 1) * sizeof(double));
    printf("T = %f\ndt_min = %f, dt_max = %f; Ndt = %d\n\n", T, dt_range[0], dt_range[Ndt-1], Ndt);

    // Creation du fichier de stockage de la complexité
    FILE *complexity_file = fopen("./data/complexity.txt", "w");
    if (!complexity_file) {
        perror("Failed to open complexity.txt");
    }

    int node_I = 1; // Nœud I


    for (int i = 0; i < Ndt; i++) {

        // Récup condition initiale --> ds les quels on va stocker les u, v au temps T
        n_m = get_intial_condition(init, u, v, N);

        dt = dt_range[i];
        int nbr_iter = (int) (T/ dt); // Nombre d'itérations de la méthode de Newmark

        int pres = 10;                  // On effectue 10 fois la méthode de Newmark pour chaque dt pour plus de précision
        int pres = 1;                   // = 1 pr lorsqu'on regarde pas la complexité temporelle
        double elapsed = 0.0;

        printf("Calcul nr = ");
        for (int j = 0; j < pres; j++){
            printf("%d / %d  \n", j+1, pres);
            clock_t start = clock();
            newmark(
                u, v,
                t,t,t,t,t,
                Ksp, Msp,
                nbr_iter, node_I,
                dt,
                n_m
            );
            clock_t end = clock(); // Fin chrono
            elapsed += ((double)(end - start) / CLOCKS_PER_SEC) / pres;
            for (int k = 0; k < N; k++) {
                u_moy[k] += u[k] / pres;
                v_moy[k] += v[k] / pres;
            }
        }
        printf("\n");


        printf("Itération %d : temps de calcul = %.8f secondes\n", i+1, elapsed);
        fprintf(complexity_file, "%f %.8f\n", dt, elapsed);

        char filename[64];
        snprintf(filename, sizeof filename, "./data/convergence/final_dt%.5g.txt", dt);
        printf("✓ dt = %.5g  -> %s (%d itérations)\n", dt, filename, nbr_iter);
        stock_final(n_m, filename, u_moy, v_moy);
    } 
    fclose(complexity_file);
    free(t);
    free(u);
    free(v);
    free(u_moy);
    free(v_moy);
}

/* ============================================================
 *  STABILITÉ : rayon spectral du schéma de Newmark
 * ============================================================*/

/* ---------- Produit y = A(h)·x  ----------------------------
 * x et y sont des vecteurs de taille 2N = [q ; v].
 * -----------------------------------------------------------*/
static void newmark_apply_amp(
        const CSRMatrix *K, const CSRMatrix *M,
        double h, double beta, double gamma,
        const double *x, double *y, int n_nodes)
{
    const int N = 2 * n_nodes;

    const double *q  = x;        /* q_n            */
    const double *v  = x + N;    /* v_n            */
          double *qn = y;        /* q_{n+1}        */
          double *vn = y + N;    /* v_{n+1}        */

    /* --------- Étape 1 : q_{n+1} --------------------------- */
    /* A_eff = M + beta h² K  */
    double *Aeff = malloc(K->nnz * sizeof(double));
    build_M_coeffK(N, (CSRMatrix *)K, (CSRMatrix *)M, Aeff,
                   beta * h * h);

    /* rhs = (M - (0.5-beta)h² K) q  +  h M v */
    double *rhs = calloc(N, sizeof(double));
    double *tmp = malloc(N * sizeof(double));

    /* tmp = M q */
    Matvec(N, M->row_ptr, M->col_idx, M->data, q, tmp);
    cblas_dcopy(N, tmp, 1, rhs, 1);

    /* rhs -= (0.5-beta)h² K q */
    if (fabs(0.5 - beta) > 1e-14) {
        Matvec(N, K->row_ptr, K->col_idx, K->data, q, tmp);
        cblas_daxpy(N, -(0.5 - beta) * h * h, tmp, 1, rhs, 1);
    }
    /* rhs += h M v */
    Matvec(N, M->row_ptr, M->col_idx, M->data, v, tmp);
    cblas_daxpy(N, h, tmp, 1, rhs, 1);

    /* q_{n+1}  =  A_eff^{-1} rhs  (CG) */
    CG(N, K->row_ptr, K->col_idx, Aeff, rhs, qn, 1e-12);

    /* --------- Étape 2 : v_{n+1} --------------------------- */
    /*  a_n      = -M^{-1} K q_n     */
    /*  a_{n+1}  = -M^{-1} K q_{n+1} */

    double *acc = malloc(N * sizeof(double));   /* réutilisé */

    Matvec(N, K->row_ptr, K->col_idx, K->data, q, tmp);
    CG(N, M->row_ptr, M->col_idx, M->data, tmp, acc, 1e-12);
    cblas_dscal(N, -1.0, acc, 1);               /* acc = a_n */

    Matvec(N, K->row_ptr, K->col_idx, K->data, qn, tmp);
    CG(N, M->row_ptr, M->col_idx, M->data, tmp, vn, 1e-12);
    cblas_dscal(N, -1.0, vn, 1);                /* vn = a_{n+1} */

    /* vn = v_n + h[ (1-gamma)a_n + gamma a_{n+1} ] */
    for (int i = 0; i < N; ++i)
        vn[i] = v[i] + h * ( (1.0 - gamma) * acc[i] + gamma * vn[i] );

    free(acc); free(tmp); free(rhs); free(Aeff);
}

/* ---------- Puissance itérée : rayon spectral --------------*/
double spectral_radius_newmark(
        const CSRMatrix *K, const CSRMatrix *M,
        int n_nodes, double h,
        double beta, double gamma,
        int    max_it, double tol)
{
    const int Ntot = 4 * n_nodes;      /* 2N avec N = 2*n_nodes        */
    double *x  = malloc(Ntot * sizeof(double));
    double *y  = malloc(Ntot * sizeof(double));

    /*  initialisation (vecteur aléatoire normalisé)  */
    srand((unsigned)time(NULL));
    for (int i = 0; i < Ntot; ++i) x[i] = (double)rand() / RAND_MAX;
    double nrm = sqrt(cblas_ddot(Ntot, x, 1, x, 1));
    cblas_dscal(Ntot, 1.0 / nrm, x, 1);

    double lambda = 0.0, lambda_old = 0.0;
    for (int it = 0; it < max_it; ++it) {

        newmark_apply_amp(K, M, h, beta, gamma, x, y, n_nodes);

        /*  rapport de normes : approximation de la valeur propre  */
        lambda = sqrt(cblas_ddot(Ntot, y, 1, y, 1));
        cblas_dscal(Ntot, 1.0 / lambda, y, 1);    /* normalise y  */

        if (fabs(lambda - lambda_old) < tol * lambda) break;
        lambda_old = lambda;

        /*  x ← y  */
        cblas_dcopy(Ntot, y, 1, x, 1);
    }

    free(x); free(y);
    return lambda;
}

/* =================================================================
 *  Balayage en pas de temps  Δt  ––  écrit spectral_radius_dt.txt
 * ================================================================= */
void stability_scan(const CSRMatrix *K, const CSRMatrix *M,
                    int n_nodes,
                    const double *dt_list, int Ndt,
                    double beta, double gamma)
{
    const char *fname = "./data/spectral_radius_dt.txt";
    FILE *fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    for (int k = 0; k < Ndt; ++k) {
        double rho = spectral_radius_newmark(
                        K, M, n_nodes,
                        dt_list[k], beta, gamma,
                        /* max_it */ 200, /* tol */ 1e-10);

        fprintf(fp, "%.15e %.15e\n", dt_list[k], rho);
        printf("dt = %-10g  |  rho = %.12g\n", dt_list[k], rho);
    }
    fclose(fp);
    puts("→ fichier « spectral_radius_dt.txt » créé.");
}

/* =================================================================
 *  Cartographie rho(β,γ) pour un pas de temps fixe
 * ================================================================= */
void stability_map_beta_gamma(const CSRMatrix *K, const CSRMatrix *M,
                              int n_nodes,
                              const double *beta_vec,  int Nbeta,
                              const double *gamma_vec, int Ngamma,
                              double dt_fixed)
{
    const char *fname = "./data/spectral_radius_bg.txt";
    FILE *fp = fopen(fname, "w");
    if (!fp) { perror(fname); return; }

    for (int i = 0; i < Nbeta; ++i)
        for (int j = 0; j < Ngamma; ++j)
        {
            double rho = spectral_radius_newmark(
                           K, M, n_nodes,
                           dt_fixed, beta_vec[i], gamma_vec[j],
                           200, 1e-10);

            /* format :  β  γ  rho  */
            fprintf(fp, "%.15e %.15e %.15e\n",
                        beta_vec[i], gamma_vec[j], rho);
        }
    fclose(fp);
    puts("→ fichier « spectral_radius_bg.txt » créé.");
}

/* ============================================================
 *  Drivers facultatifs pour tester rapidement depuis main()
 * ============================================================*/
void run_stability_vs_dt(const CSRMatrix *K, const CSRMatrix *M,
                         int n_nodes)
{
    static const double dt_range[] = {0.0001, 0.001, 0.002, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0,
                            2.0, 5.0, 10.0, 15.0, 20.0, 50.0, 100.0};
    int Ndt = (int)(sizeof dt_range / sizeof dt_range[0]);

    stability_scan(K, M, n_nodes,
                   dt_range, Ndt,
                   0.25, 0.50);          /* β = ¼  γ = ½ */
}

void run_stability_vs_bg(const CSRMatrix *K, const CSRMatrix *M,
                         int n_nodes)
{
    static const double beta_tab[]  = {0.05,0.10,0.15,0.20,0.25,0.30,0.35,
                                       0.40,0.45,0.50};
    static const double gamma_tab[] = {0.1,0.20,0.30,0.40,0.50,0.60,0.7,0.8,0.9,1};

    stability_map_beta_gamma(K, M, n_nodes,
                             beta_tab,  (int)(sizeof beta_tab  / sizeof beta_tab[0]),
                             gamma_tab, (int)(sizeof gamma_tab / sizeof gamma_tab[0]),
                             0.10);                 /* h fixé = 0.10 */
}