#include "devoir_2.h"
#include "utils.h"
#include "model.h"
#include "utils_gmsh.h"
#include "gmshc.h"
#include "../gmsh-sdk/include/gmshc.h"
#include <math.h>
#include <cblas.h>
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
        printf("Usage: \n./deformation <model> <mesh_size_ratio>\n");
        printf("model: one of the model implemented in models/\n");
        printf("mesh_size_ratio: mesh size factor\n");
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
    CG(Ksp->n, Ksp->row_ptr, Ksp->col_idx, Ksp->data, rhs, sol, eps);    
    //display_sol(model, sol);

    // Début modif
    // Devoir 3

    // Récup ux_i uy_i vx_i vy_i de <initial.txt>
    int n = Kbd->n;

    double *uxi = (double *)malloc( n * sizeof(double));
    double *uyi = (double *)malloc( n * sizeof(double));
    double *uvxi = (double *)malloc( n * sizeof(double));
    double *uvyi = (double *)malloc( n * sizeof(double));
    n = get_intial_condition(argv[5], uxi, uyi, uvxi, uvyi);
    printf("n  = %d\n", n);


    // Newmark computation
    double *t = (double *)malloc( n * sizeof(double));
    //...


    // Stockage solution dans <final.txt> : le déplacement et la vitesse au temps T, au même format que le fichier <initial.txt>
    stock_final(n, argv[6], uxi, uyi, uvxi, uvyi);
    

    // Stockage dans <time.txt> le déplacement et la vitesse d’un nœud I à chaque itération temporelle
    stock_time(n, argv[7], t, uxi, uyi, uvxi, uvyi);
    


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
