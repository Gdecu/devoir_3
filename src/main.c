#include "devoir_2.h"
#include "utils.h"
#include "model.h"
#include "utils_gmsh.h"
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

    // Début modif
    // Devoir 3

    // Récup ux_i uy_i vx_i vy_i de <initial.txt>
    int nnz = Ksp->nnz;
    int n = Ksp->n;
    printf("nnz  = %d\n", nnz);
    printf("n  = %d\n", n);

    // Récup condition initiale --> ds les quels on va stocker les u, v au temps T
    double *u = (double *)malloc( n * sizeof(double));
    double *v = (double *)malloc( n * sizeof(double));
    n = get_intial_condition(argv[5], u, v, n);
    printf("ux[0]  = %15le, uy[0] = %15le\n", u[0], u[1]);
    printf("vx[0]  = %15le, vy[0] = %15le\n", v[0], v[1]);
    printf("n  = %d\n", n);
    printf("\n");

    int T = atoi(argv[3]); // Temps total
    double dt = atof(argv[4]); // Pas de temps
    printf("T  = %d\n", T);
    printf("dt  = %le\n", dt);
    int node_I = atoi(argv[8]); // Nœud I

    if (T <= 0 || dt <= 0) {
        printf("Erreur : T et dt doivent être des valeurs positives.\n");
        return -1;
    }

    // On vas stocker les uxI uyI vxI vyI du nœud I à chaque itération temporelle
    int N = (int) round(T / dt);
    double *uxI = (double *)malloc( N * sizeof(double));
    double *uyI = (double *)malloc( N * sizeof(double));
    double *vxI = (double *)malloc( N * sizeof(double));
    double *vyI = (double *)malloc( N * sizeof(double));
    double *t = (double *)malloc( N * sizeof(double));


    // Newmark computation
    newmark(
        u, v,
        uxI, uyI, vxI, vyI, t,
        Ksp, Msp,
        T, node_I,
        dt,
        n
    );


    // Stockage solution dans <final.txt> : le déplacement et la vitesse au temps T, au même format que le fichier <initial.txt>
    stock_final(n, argv[6], u, v);
    

    // Stockage dans <time.txt> le déplacement et la vitesse d’un nœud I à chaque itération temporelle
    stock_time(N, argv[7], t, uxI, uyI, vxI, vyI);



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
