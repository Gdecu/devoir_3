#include "devoir_2.h"
#include "utils.h"
#include "model.h"
#include "utils_gmsh.h"
#include "analyse.h"
#include "gmshc.h"
#include "../gmsh-sdk/include/gmshc.h"
#include <math.h>
#include <gmshc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define VERBOSE 1
#define PRECISION 10

void display_sol(FE_Model *model, double *sol) {
    int ierr, n_views, *views;
    double *bounds;
    add_gmsh_views(&views, &n_views, &bounds);

    double *data_forces = malloc(6 * model->n_bd_edge * sizeof(double));
    visualize_disp(model, sol, views[1], 0, &bounds[2]);
    visualize_stress(model, sol, views, 1, 0, data_forces, bounds);
    visualize_bd_forces(model, data_forces, views[0], 1, &bounds[0]);

    create_tensor_aliases(views);
    set_view_options(n_views, views, bounds);
    gmshFltkRun(&ierr);
    gmshFltkFinalize(&ierr);
    free(data_forces);
}

void display_info(FE_Model *model, int step, struct timespec ts[4]) {

    char *m_str[3] = {"Plane stress", "Plane strain", "Axisymmetric"};
    char *r_str[4] = {"No", "X", "Y", "RCMK"};

    if (step == 1) {
        printf(
            "\n===========  Linear elasticity simulation - FEM  ===========\n\n"
        );
        printf("%30s = %s\n", "Model", model->model_name);
        printf("%30s = %s\n", "Model type", m_str[model->m_type]);
        printf("%30s = %.3e\n", "Young's Modulus E", model->E);
        printf("%30s = %.3e\n", "Poisson ratio nu", model->nu);
        printf("%30s = %.3e\n\n", "Density rho", model->rho);
    } else if (step == 2) {
        char *e_str = (model->e_type == TRI) ? "Triangle" : "Quadrilateral";
        printf("%30s = %s\n", "Element type", e_str);
        printf("%30s = %zu\n", "Number of elements", model->n_elem);
        printf("%30s = %zu\n", "Number of nodes", model->n_node);
        printf("%30s = %s\n", "Renumbering", r_str[model->renum]);
        printf("%30s = %zu\n\n", "Matrix bandwidth", 2 * model->node_band + 1);
    }
}


int main(int argc, char *argv[]) {

    int ierr;
    double mesh_size_ratio;
    if ((argc < 3) || (sscanf(argv[2], "%lf", &mesh_size_ratio)) != 1) {
        printf("Usage: \n./deformation <model> <lc> <T> <dt> <initial.txt> <final.txt> <time.txt> <I>\n");
        printf("model: one of the model implemented in models/\n");
        return -1;
    }

    // Simulation parameters
    const ElementType e_type = TRI;
    const Renumbering renum = RENUM_NO;  // let gmsh do the RCMK renumbering

    FE_Model *model = create_FE_Model(argv[1], e_type, renum);
    display_info(model, 1, NULL);

    gmshInitialize(argc, argv, 0, 0, &ierr);
    gmshOptionSetNumber("General.Verbosity", 2, &ierr);
    model->mesh_model(mesh_size_ratio, e_type);

    load_mesh(model);
    renumber_nodes(model);
    display_info(model, 2, NULL);
    assemble_system(model);
    double *rhs = (double *)calloc(2 * model->n_node, sizeof(double));
    double *sol = (double *)calloc(2 * model->n_node, sizeof(double));
    add_bulk_source(model, rhs);
    enforce_bd_conditions(model, rhs);

    SymBandMatrix *Kbd = model->K;
    SymBandMatrix *Mbd = model->M;
    CSRMatrix *Ksp = band_to_sym_csr(Kbd);
    CSRMatrix *Msp = band_to_sym_csr(Mbd);
    double eps = 1e-8;
    //CG(Ksp->n, Ksp->row_ptr, Ksp->col_idx, Ksp->data, rhs, sol, eps);    
    //display_sol(model, sol);

    //
    // Début de l'implémentation du devoir 3 : la méthode de Newmark
    //

    printf("==========================================================================================================================\n");
    printf("Début de l'implémentation du devoir 3 : la méthode de Newmark\n");
    printf("==========================================================================================================================\n\n");
    // Récup ux_i uy_i vx_i vy_i de <initial.txt>
    int nnz = Ksp->nnz;
    int n = Ksp->n;
    printf("nnz  = %d\n", nnz);
    printf("n  = %d\n\n", n);

    // Récup condition initiale --> ds les quels on va stocker les u, v au temps T
    double *u = (double *)malloc( n * sizeof(double));
    double *v = (double *)malloc( n * sizeof(double));
    printf("Recup condition initiale\n\n");
    n = get_intial_condition(argv[5], u, v, n);
    //printf("ux[0]  = %15le, uy[0] = %15le\n", u[0], u[1]);
    //printf("vx[0]  = %15le, vy[0] = %15le\n", v[0], v[1]);
    printf("n  = %d\n", n);
    printf("\n");

    int T = atoi(argv[3]); // Temps total
    double dt = atof(argv[4]); // Pas de temps
    int nbr_iter = (int) (T/ dt); // Nombre d'itérations de la méthode de Newmark
    printf("T  = %d\n", T);
    printf("dt  = %le\n", dt);
    printf("nbr_iter  = %d\n\n", nbr_iter);
    int node_I = atoi(argv[8]); // Nœud I

    if (T <= 0 || dt <= 0) {
        printf("Erreur : T et dt doivent être des valeurs positives.\n");
        return -1;
    }

    // On vas stocker les uxI uyI vxI vyI du nœud I à chaque itération temporelle
    double *uxI = (double *)malloc( (nbr_iter + 1) * sizeof(double));
    double *uyI = (double *)malloc( (nbr_iter + 1) * sizeof(double));
    double *vxI = (double *)malloc( (nbr_iter + 1) * sizeof(double));
    double *vyI = (double *)malloc( (nbr_iter + 1) * sizeof(double));
    double *t = (double *)malloc( (nbr_iter + 1) * sizeof(double));

    uxI[0] = u[2 * node_I];
    uyI[0] = u[2 * node_I + 1];
    vxI[0] = v[2 * node_I];
    vyI[0] = v[2 * node_I + 1];
    t[0]   = 0.0;


    // Newmark computation
    printf("Début de la méthode de Newmark...\n");
    newmark(
        u, v,
        uxI, uyI, vxI, vyI, t,
        Ksp, Msp,
        nbr_iter, node_I,
        dt,
        n
    );
    printf("Fin de la méthode de Newmark\n\n");


    // Stockage solution dans <final.txt> : le déplacement et la vitesse au temps T, au même format que le fichier <initial.txt>
    printf("Stockage solution dans <final.txt> : le déplacement et la vitesse au temps T, au même format que le fichier <initial.txt>\n");
    stock_final(n, argv[6], u, v);


    // Stockage dans <time.txt> le déplacement et la vitesse d’un nœud I à chaque itération temporelle
    printf("Stockage dans <time.txt> le déplacement et la vitesse d’un nœud I à chaque itération temporelle\n");
    stock_time(nbr_iter+1, argv[7], t, uxI, uyI, vxI, vyI);

    printf("==========================================================================================================================\n");
    printf("Fin de l'implémentation du devoir 3 : la méthode de Newmark\n");
    printf("==========================================================================================================================\n\n");

    //
    // Fin de l'implémentation du devoir 3 : la méthode de Newmark
    //

    //
    // Analyses
    //
    printf("==========================================================================================================================\n");
    printf("Début des analyses\n");
    printf("==========================================================================================================================\n\n");

    // Récup les coords du modèle (utile pour l'animation)
    printf("Récup les coords du modèle (utile pour l'animation)\n\n");
    char *filename = "./data/coords.txt";
    get_coords(model->coords, model->n_node, filename);

    // Verification que l'énergie cinétique et potentielle sont conservées
    printf("Verification que l'énergie cinétique et potentielle sont conservées\n");
    // Animation
    printf("Récupreration de %d etats final pour la simulation\n\n", nbr_iter);
    T = 100;
    //anim_enrgy(argv[5], Ksp, Msp, u, v, n, T, dt, nbr_iter);
    printf("\n\n");

    // Test ordre de convergence & complexité temporelle  -- on fait varier dt a T cst
    printf("Analyse ordre de convergence\n");
    printf("Analyse complexité temporelle\n");

    // Create sol_ref
    double dt_range_ref[1] = {0.0001};
    int Ndt_ref = 1;
    T = 200;
    //convergence_complexity(Ksp, Msp, 2*n, T, argv[5], dt_range_ref, Ndt_ref);

    // Create solution for a range of dt
    int Ndt = 22;
    double dt_range[22] = { 0.001, 0.002, 0.0025, 0.005, 0.0075, 0.01, 0.015,
                            0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0,
                            2.0, 5.0, 10.0, 15.0, 20.0, 50.0, 100.0};
    T = 200;
    //convergence_complexity(Ksp, Msp, 2*n, T, argv[5], dt_range, Ndt);
    printf("\n\n");

    // Test stabilité du schéma de Newmark
    printf("Test stabilité du schéma de Newmark : rayon spectral\n");
    run_stability_vs_dt(Ksp, Msp, n);
    run_stability_vs_bg(Ksp, Msp, n);

    printf("==========================================================================================================================\n");
    printf("Fin des analyses\n");
    printf("==========================================================================================================================\n\n");


    free(u);
    free(v);

    free(uxI);
    free(uyI);
    free(vxI);
    free(vyI);
    free(t);
    // fin modif


    
    // Free stuff
    free_csr(Ksp);
    free_csr(Msp);
    gmshFinalize(&ierr);
    free(sol);
    free(rhs);
    free_FE_Model(model);
    return 0;
}
