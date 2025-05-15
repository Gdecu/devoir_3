#ifndef ANALYSE_H
#define ANALYSE_H

#include "model.h"

// Stuff for animation
void get_nbrIter_finalSol(char *initial_conditions, CSRMatrix *Ksp, CSRMatrix *Msp, double *u, double *v, int n, int T, double dt, int nbr_iter);
void display_anim();

#endif