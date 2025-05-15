#include "devoir_2.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "model.h"
#include "analyse.h"

void newmark_anim(
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
    
    for (int k = 0; k < nbr_iter; k++) {

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

        // Stockage solution dans <final.txt> : le déplacement et la vitesse au temps T, au même format que le fichier <initial.txt>
        char filename[40];
        sprintf(filename, "./data/anim/final_%d.txt", k+1);
        stock_final(n_nodes, filename, u, v);
    }

    free(p); free(rhs); free(q_new); free(Kq); free(q_moy); free(Beff);
    free(Aeff);
}

void get_nbrIter_finalSol(char *initial_conditions, CSRMatrix *Ksp, CSRMatrix *Msp, double *u, double *v, int n, int T, double dt, int nbr_iter) {

    // On stocke l'état intiale
    get_intial_condition(initial_conditions, u, v, 2*n);
    //printf("ux[0]  = %15le, uy[0] = %15le\n", u[0], u[1]);
    //printf("vx[0]  = %15le, vy[0] = %15le\n", v[0], v[1]);
    stock_final(n, "./data/anim/final_0.txt", u, v);

    int nnz = Ksp->nnz;
    printf("K->nnz  = %d\n", nnz);
    printf("n  = %d\n\n", n);

    printf("T  = %d\n", T);
    printf("dt  = %le\n", dt);
    printf("nbr_iter  = %d\n\n", nbr_iter);

    newmark_anim(u, v,
        initial_conditions,
        Ksp, Msp, 
        nbr_iter, dt , n);
}

void display_anim(){

}