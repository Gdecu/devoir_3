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

int get_intial_condition(const char *filename, double *uxi, double *uyi, double *uvxi, double *uvyi);
void stock_final(int n, const char *filename, double *uxi, double *uyi, double *uvxi, double *uvyi);
void stock_time(int n, const char *filename, double *t, double *uxi, double *uyi, double *uvxi, double *uvyi);

#endif
