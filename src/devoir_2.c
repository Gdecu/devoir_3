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
    double *uxi, double *uyi, double *vxi, double *vyi,
    double *uxI, double *uyI, double *vxI, double *vyI, double *t,
    double *Kval, double *M,
    double T_final, int node_I,
    double dt,
    int n
) {
    int T = (int)(T_final / dt);
    double beta = 0.25;
    double gamma = 0.5;

    extern int *rows_idx, *cols;
    double *A_eff = (double *)malloc(rows_idx[n] * sizeof(double));

    build_effective_matrix(n, rows_idx, cols, Kval, M, A_eff, beta, dt);

    for (int i = 0; i < T; i++) {
        // Itération Newmark pour chaque direction
        newmark_iter(uxi, M, Kval, A_eff, rows_idx, cols, dt, beta, gamma, n);
        newmark_iter(uyi, M, Kval, A_eff, rows_idx, cols, dt, beta, gamma, n);
        newmark_iter(vxi, M, Kval, A_eff, rows_idx, cols, dt, beta, gamma, n);
        newmark_iter(vyi, M, Kval, A_eff, rows_idx, cols, dt, beta, gamma, n);

        // Stockage du noeud d'intérêt
        uxI[i] = uxi[node_I];
        uyI[i] = uyi[node_I];
        vxI[i] = vxi[node_I];
        vyI[i] = vyi[node_I];
        t[i] = i * dt;
    }

    free(A_eff);
}


void newmark_iter(
    double *q, // uxi / uyi / vxi / vyi
    double *M,
    double *K,
    double *A_eff, // matrice M + β h^2 K
    int *rows_idx,
    int *cols,
    double h,
    double beta,
    double gamma,
    int n
) {
    double *p = (double *)malloc(n * sizeof(double));
    double *rhs = (double *)malloc(n * sizeof(double));
    double *qn1 = (double *)malloc(n * sizeof(double));
    double *Kq = (double *)malloc(n * sizeof(double));
    double *tmp = (double *)malloc(n * sizeof(double));

    // p_n = M * q̇_n
    for (int i = 0; i < n; i++) {
        p[i] = M[i] * q[i]; // q = q̇ ici
    }

    // K * q_n
    Matvec(n, rows_idx, cols, K, q, Kq);

    // RHS = (M - h²/2 (1 - 2β) K) * q + h * p
    for (int i = 0; i < n; i++) {
        rhs[i] = (M[i] - h * h / 2 * (1 - 2 * beta) * Kq[i]) * q[i] + h * p[i];
    }

    // Résolution A_eff * q_{n+1} = rhs
    CG(n, rows_idx, cols, A_eff, rhs, qn1, 1e-10);

    // tmp = (1 - γ) * q + γ * qn1
    for (int i = 0; i < n; i++) {
        tmp[i] = (1 - gamma) * q[i] + gamma * qn1[i];
    }

    Matvec(n, rows_idx, cols, K, tmp, Kq);

    // p_{n+1} = p - h * K * tmp
    for (int i = 0; i < n; i++) {
        p[i] -= h * Kq[i];
    }

    // q_{n+1} = qn1, q̇_{n+1} = p / M
    for (int i = 0; i < n; i++) {
        q[i] = qn1[i];
        q[i] = p[i] / M[i]; // q est maintenant q̇
    }

    free(p); free(rhs); free(qn1); free(Kq); free(tmp);
}