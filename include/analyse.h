#ifndef ANALYSE_H
#define ANALYSE_H

#include "model.h"

#define MAXNAME 256
#define FALSE 0 
#define TRUE  1

// Stuff for animation
void anim_enrgy(char *initial_conditions, CSRMatrix *Ksp, CSRMatrix *Msp, double *u, double *v, int n, int T, double dt, int nbr_iter);
void get_coords(double *coord, int n_nodes, char *filename);
void convergence_complexity(CSRMatrix *Ksp, CSRMatrix *Msp, int N, double T, char *init, double *dt_range, int Ndt);
#endif