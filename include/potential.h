#ifndef ESF_POTENTIAL_H
#define ESF_POTENTIAL_H
#include "tube.h"
typedef struct {
    int first_row, last_row, nx, ny, n, nz, neumann_x, periodic_x, *row, *col, *diag;
    double *a, *lu, *rhs, *phi, *work;
    double residual;
    int iterations;
} Potential;
int potential_create(Potential *p, int nx, int ny);
int potential_create_periodic(Potential *p, int nx, int ny);
void potential_free(Potential *p);
void potential_zero(Potential *p);
void potential_element(Potential *p, int i, int j, const double K[4], const double C[2]);
void potential_boundary(Potential *p);
int potential_solve(Potential *p, double rtol);
void potential_coefficients(const TubeGrid *g, const double *rho, int first, int count, int hall,
                            double *coefficients);
void potential_assemble(Potential *p, const double *coefficients);
#endif
