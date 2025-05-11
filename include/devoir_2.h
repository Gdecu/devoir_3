#ifndef DEVOIR_2_H
#define DEVOIR_2_H

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

int get_intial_condition(const char *filename, double *uxi, double *uyi, double *vxi, double *vyi);
void stock_final(int n, const char *filename, double *uxi, double *uyi, double *vxi, double *vyi);
void stock_time(int T, const char *filename, double *t, double *uxi, double *uyi, double *vxi, double *vyi);

void newmark(
    double *uxi, double *uyi, double *vxi, double *vyi,
    double *uxI, double *uyI, double *vxI, double *vyI, double *t,
    double *K, double *M,
    double T, int node_I,
    double dt,
    int n
);

void newmark_iter(
    double *uxi, double *uyi, double *vxi, double *vyi,
    double *M, double *K,
    int n
);
#endif
