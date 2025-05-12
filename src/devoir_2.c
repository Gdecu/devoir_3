#include "devoir_2.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// SYM = 1 --> Vous ne stockez que la partie inferieure de la matrice
// SYM = 0 --> Vous stockez toutes les entrées non nulles de la matrice
#define SYM 1 
int csr_sym() { return SYM; }

/**
 * Si vous avez installé la librairie BLAS, vous pouvez l'ajouter au Makefile
 * et supprimer les implémentations de cblas_* ci-dessous.
 */

// ----------------------- BLAS start -----------------------

void cblas_dscal(int n, double alpha, double *x, int incx) {
    for (int i = 0; i < n; i++) {
        x[i * incx] *= alpha;
    }
}
void cblas_dcopy(int n, const double *x, int incx, double *y, int incy) {
    for (int i = 0; i < n; i++) {
        y[i * incy] = x[i * incx];
    }
}
double cblas_ddot(int n, const double *x, int incx, const double *y, int incy) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += x[i * incx] * y[i * incy];
    }
    return result;
}
void cblas_daxpy(
    int n, double alpha, const double *x, int incx, double *y, int incy
) {
    for (int i = 0; i < n; i++) {
        y[i * incy] += alpha * x[i * incx];
    }
}
double cblas_dnrm2(int n, const double *x, int incx) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += x[i * incx] * x[i * incx];
    }
    return sqrt(result);
}

// ------------------------ BLAS end ------------------------

void Matvec(
    int n,
    const int *rows_idx,
    const int *cols,
    const double *A,
    const double *v,
    double *Av
) {

    cblas_dscal(n, 0, Av, 1);

    for (int i = 0; i < n; i++) {
#if SYM
        for (int j = rows_idx[i]; j < rows_idx[i + 1] - 1; j++) {
            Av[i] += A[j] * v[cols[j]];
            Av[cols[j]] += A[j] * v[i];
        }
        Av[i] += A[rows_idx[i + 1] - 1] * v[i];
#else
        for (int j = rows_idx[i]; j < rows_idx[i + 1]; j++) {
            Av[i] += A[j] * v[cols[j]];
        }
#endif
    }
}

int CG(
    int n,
    const int *rows_idx,
    const int *cols,
    const double *A,
    const double *b,
    double *x,
    double eps
) {
    int idx = 0;
    double alpha, beta;
    double *p = (double *)malloc(n * sizeof(double));
    double *Ap = (double *)malloc(n * sizeof(double));
    double *r = (double *)malloc(n * sizeof(double));

    cblas_dscal(n, 0, x, 1);
    cblas_dcopy(n, b, 1, r, 1);
    cblas_dcopy(n, r, 1, p, 1);
    double r0 = cblas_dnrm2(n, b, 1);
    double r_norm2 = cblas_ddot(n, r, 1, r, 1);
    // printf("r0 = %9.3le\n", r0);

    while (sqrt(r_norm2) / r0 > eps) {
        // if (idx %100 == 0) printf("idx : %4d\n", idx);
        Matvec(n, rows_idx, cols, A, p, Ap);
        alpha = r_norm2 / cblas_ddot(n, p, 1, Ap, 1);
        cblas_daxpy(n, alpha, p, 1, x, 1);
        cblas_daxpy(n, -alpha, Ap, 1, r, 1);
        beta = 1 / r_norm2;
        r_norm2 = cblas_ddot(n, r, 1, r, 1);
        beta *= r_norm2;
        cblas_dscal(n, beta, p, 1);
        cblas_daxpy(n, 1, r, 1, p, 1);
        idx++;
        // printf("it : %3d  -> res = %9.3le\n", idx, sqrt(r_norm2) / r0);
    }

    // printf("cg it : %d\n", idx);
    free(p);
    free(Ap);
    free(r);
    return idx;
}

//
// debut modif
//

// Récup ux_i uy_i vx_i vy_i de <initial.txt>
int get_intial_condition(
    const char *filename,
    double *uxi,
    double *uyi,
    double *vxi,
    double *vyi
) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error opening file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    int i = 0;
    while (fscanf(fp, "%lf %lf %lf %lf", &uxi[i], &uyi[i], &vxi[i], &vyi[i]) == 4) {
        i++;
    }

    fclose(fp);
    return i;
}


// Stockage solution dans <final.txt> : le déplacement et la vitesse au temps T, au même format que le fichier <initial.txt>
void stock_final(int n,
    const char *filename,
    double *uxi,
    double *uyi,
    double *vxi,
    double *vyi
) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error opening file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; i++) {
        fprintf(fp, "%lf %lf %lf %lf\n", uxi[i], uyi[i], vxi[i], vyi[i]);
    }

    fclose(fp);
}


// Stockage dans <time.txt> le déplacement et la vitesse d’un nœud I à chaque itération temporelle
void stock_time(
    int T,
    const char *filename,
    double *t,
    double *uxi,
    double *uyi,
    double *vxi,
    double *vyi
) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error opening file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < T; i++) {
        fprintf(fp, "%lf %lf %lf %lf %lf\n",t[i], uxi[i], uyi[i], vxi[i], vyi[i]);
    }

    fclose(fp);
}

void build_effective_matrix(
    int n,
    const int *rows_idx,
    const int *cols,
    const double *Kval,
    const double *M,
    double *A_eff,
    double beta,
    double h
) {
    double coeff = beta * h * h;

    for (int i = 0; i < n; i++) {
        for (int j = rows_idx[i]; j < rows_idx[i + 1]; j++) {
            int col = cols[j];
            A_eff[j] = coeff * Kval[j];
            if (i == col) {
                A_eff[j] += M[i]; 
            }
        }
    }
}



void newmark(
    double *uxi, double *uyi,      // initialisations u
    double *vxi, double *vyi,      // initialisations v
    double *uxI, double *uyI,      // sorties du nœud I
    double *vxI, double *vyI,
    double *t,                      // temps
    double *Kval, double *Mdiag,    // CSR K et vecteur M
    double Tfinal, int node_I,
    double dt, int n_nodes
) {
    int N    = 2 * n_nodes;                // DOF
    int T    = (int)round(Tfinal / dt);
    double β = 0.25, γ = 0.5;

    extern int *rows_idx, *cols;           // définis par ta lecture CSR
    double *Aeff = malloc(rows_idx[N] * sizeof(double));
    build_effective_matrix(N, rows_idx, cols, Kval, Mdiag, Aeff, β, dt);

    // empile uxi/uyi et vxi/vyi dans q et v
    double *q = malloc(N * sizeof(double));
    double *v = malloc(N * sizeof(double));
    for (int i = 0; i < n_nodes; i++) {
        q[2*i]   = uxi[i];
        q[2*i+1] = uyi[i];
        v[2*i]   = vxi[i];
        v[2*i+1] = vyi[i];
    }

    for (int k = 0; k < T; k++) {
        newmark_iter(N, rows_idx, cols,
                     Kval, Aeff, Mdiag,
                     q, v, dt, β, γ);

        // extrait la DOF du nœud_I
        uxI[k] = q[2*node_I];
        uyI[k] = q[2*node_I+1];
        vxI[k] = v[2*node_I];
        vyI[k] = v[2*node_I+1];
        t[k]   = k * dt;
    }

    free(Aeff);
    free(q);
    free(v);
}



/**
 * Mise à jour q ← q_{n+1} et v ← v_{n+1} en une seule passe,
 * pour le système non amorti.
 */
void newmark_iter(
    int N,                    // nombre total de DOF (2·#nœuds)
    const int *rows_idx,
    const int *cols,
    const double *Kval,       // K en CSR
    const double *Aeff,       // M + βh²K en CSR
    const double *Mdiag,      // masse diagonale
    double *q,                // [ux,uy,…], taille N
    double *v,                // vitesses, taille N
    double h,
    double beta,
    double gamma
) {
    double *p    = malloc(N * sizeof(double));
    double *rhs  = malloc(N * sizeof(double));
    double *qnp1 = malloc(N * sizeof(double));
    double *Kq   = malloc(N * sizeof(double));
    double *tmp  = malloc(N * sizeof(double));

    // 1) p_n = M·v_n
    for (int i = 0; i < N; i++) p[i] = Mdiag[i] * v[i];

    // 2) Kq = K·q_n
    Matvec(N, rows_idx, cols, Kval, q, Kq);

    // 3) construction du second membre
    //    rhs[i] = (Mdiag[i] - ½h²(1-2β)·Kii) q[i] + h·p[i],
    //    mais on a Kq[i] = (K·q_n)[i], donc :
    for (int i = 0; i < N; i++) {
        rhs[i] = (Mdiag[i]
                  - 0.5 * h*h * (1.0 - 2.0*beta) * Kq[i] / q[i]
                 ) * q[i]
               + h * p[i];
    }

    // 4) résolution Aeff · q_{n+1} = rhs
    CG(N, rows_idx, cols, Aeff, rhs, qnp1, 1e-10);

    // 5) p_{n+1} = p_n - h·K·[ (1-γ)q_n + γ q_{n+1} ]
    for (int i = 0; i < N; i++)
        tmp[i] = (1.0 - gamma) * q[i] + gamma * qnp1[i];
    Matvec(N, rows_idx, cols, Kval, tmp, Kq);
    for (int i = 0; i < N; i++)
        p[i] -= h * Kq[i];

    // 6) mise à jour finale : q ← q_{n+1}, v ← p / Mdiag
    for (int i = 0; i < N; i++) {
        q[i] = qnp1[i];
        v[i] = p[i] / Mdiag[i];
    }

    free(p); free(rhs); free(qnp1); free(Kq); free(tmp);
}
