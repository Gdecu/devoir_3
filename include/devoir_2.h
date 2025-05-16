#ifndef DEVOIR_2_H
#define DEVOIR_2_H

#include "model.h"

void Matvec(
    int n,
    const int *rows_idx,
    const int *cols,
    const double *A,
    const double *v,
    double *Av
);

int CG(
    int n,
    const int *rows_idx,
    const int *cols,
    const double *A,
    const double *b,
    double *x,
    double eps
);

int csr_sym();

int get_intial_condition(const char *filename, double *u, double *v, int n);
void stock_final(int n, const char *filename, double *u, double *v);
void stock_time(int T, const char *filename, double *t, double *uxi, double *uyi, double *vxi, double *vyi);

void newmark(
    double *u, double *v,          // initialisations u et v
    double *uxI, double *uyI,      // sorties du nœud I
    double *vxI, double *vyI,
    double *t,                      // temps
    CSRMatrix *K, CSRMatrix *M,    // CSR K et vecteur M
    double Tfinal, int node_I,
    double dt, int n_nodes
);

void cblas_daxpy(int n, double alpha, const double *x, int incx, double *y, int incy);
double cblas_ddot(int n, const double *x, int incx, const double *y, int incy);
double cblas_dnrm2(int n, const double *x, int incx);
void cblas_dscal(int n, double alpha, double *x, int incx);
void cblas_dcopy(int n, const double *x, int incx, double *y, int incy);
void build_M_coeffK(
    int n,
    CSRMatrix *K,
    CSRMatrix *M,
    double *A_eff,
    double coeff_K
);
#endif
