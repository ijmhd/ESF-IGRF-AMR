#include "solver.h"
#include <petscksp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Check distributed PETSc assembly/solution against the uniform reference problem.
 */
int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double worst = 0;
    for (int mode = 0; mode < 3; mode++) {
        Solver3 s = {.rank = rank, .size = size};
        s.c.grid.nx = 32;
        s.c.grid.ny = 24;
        s.c.rtol = 1e-11;
        s.c.periodic_x = mode == 1;
        Potential ref;
        if (mode == 1) {
            potential_create_periodic(&s.p, 32, 24);
            potential_create_periodic(&ref, 32, 24);
        } else {
            potential_create(&s.p, 32, 24);
            potential_create(&ref, 32, 24);
        }
        ref.neumann_x = s.p.neumann_x = mode == 2;
        s.coeff = calloc(32 * 24 * 6, sizeof(double));
        s.phi = calloc(s.p.n + 1, sizeof(double));
        for (int j = 0; j < 24; j++)
            for (int i = 0; i < 32; i++) {
                double x = (i + .5) / 32., y = (j + .5) / 24., pi = 3.141592653589793;
                double sigma = .01 + exp(-20 * (y - .5) * (y - .5)),
                       hall = .35 * sigma * sin(2 * pi * x);
                double *a = s.coeff + 6 * (j * 32 + i);
                a[0] = 1.4 * sigma;
                a[3] = .7 * sigma;
                a[1] = .2 * sigma + hall;
                a[2] = .2 * sigma - hall;
                a[4] = cos(2 * pi * x) * sin(pi * y);
                a[5] = .3 * sin(2 * pi * x) * cos(pi * y);
            }
        potential_assemble(&ref, s.coeff);
        if (potential_solve(&ref, 1e-12) < 0)
            fail3(&s, "reference solve failed");
        if (solver3_petsc_create(&s) || solver3_petsc_solve(&s))
            fail3(&s, "parallel solve failed");
        double *matrix = calloc(s.p.nz, sizeof(double)), *load = calloc(s.p.n, sizeof(double));
        MPI_Allreduce(s.p.a, matrix, s.p.nz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(s.p.rhs, load, s.p.n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int k = 0; k < s.p.nz; k++)
            if (matrix[k] != ref.a[k])
                fail3(&s, "distributed matrix differs");
        for (int k = 0; k < s.p.n; k++)
            if (load[k] != ref.rhs[k])
                fail3(&s, "distributed load differs");
        double scale = 1, err = 0;
        for (int k = 0; k < s.p.n; k++)
            scale = fmax(scale, fabs(ref.phi[k]));
        for (int k = 0; k < s.p.n; k++)
            err = fmax(err, fabs(s.phi[k] - ref.phi[k]) / scale);
        if (err > 2e-8 || s.p.residual > 1.05 * s.c.rtol)
            fail3(&s, "parallel solution differs");
        worst = fmax(worst, err);
        if (!rank)
            printf(
                "mode=%d matrix_exact=1 solution_relative_error=%.9g true_residual=%.9g iterations=%d\n",
                mode, err, s.p.residual, s.p.iterations);
        solver3_petsc_destroy(&s);
        potential_free(&s.p);
        potential_free(&ref);
        free(s.coeff);
        free(s.phi);
        free(matrix);
        free(load);
    }
    if (!rank)
        printf("PETSc distributed matrix/solution checks passed ranks=%d worst=%.9g\n", size,
               worst);
    PetscFinalize();
    return 0;
}
