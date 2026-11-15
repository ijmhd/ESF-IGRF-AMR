/*
 * Small/reference sparse transverse potential solver with bilinear elements and BiCGSTAB/ILU.
 * The production AMR path uses the separate distributed PETSc implementation.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "potential.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
/**
 * Classify a potential node or apply its prescribed boundary role; inspect the explicit
 * regional/periodic branches below.
 */
static int boundary(const Potential *p, int i) {
    return i / (p->nx + 1) == 0 || i / (p->nx + 1) == p->ny ||
           (!p->neumann_x && !p->periodic_x && (i % (p->nx + 1) == 0 || i % (p->nx + 1) == p->nx));
}
/**
 * Map a periodic or physical transverse node to its canonical independent index.
 */
static int canonical(const Potential *p, int i) {
    return p->periodic_x && i % (p->nx + 1) == p->nx ? i - p->nx : i;
}
/**
 * Allocate the reference sparse potential graph and solver work arrays for the requested
 * periodic policy.
 */
static int create(Potential *p, int nx, int ny, int periodic) {
    memset(p, 0, sizeof(*p));
    p->nx = nx;
    p->ny = ny;
    p->n = (nx + 1) * (ny + 1);
    p->first_row = 0;
    p->last_row = p->n;
    p->periodic_x = periodic;
    p->row = calloc(p->n + 1, sizeof(int));
    p->col = calloc(p->n * 9, sizeof(int));
    p->diag = calloc(p->n, sizeof(int));
    p->a = calloc(p->n * 9, sizeof(double));
    p->lu = calloc(p->n * 9, sizeof(double));
    p->rhs = calloc(p->n, sizeof(double));
    p->phi = calloc(p->n, sizeof(double));
    p->work = calloc((size_t)p->n * 10, sizeof(double));
    if (!p->row || !p->col || !p->diag || !p->a || !p->lu || !p->rhs || !p->phi || !p->work)
        return -1;
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++) {
            int r = j * (nx + 1) + i;
            p->row[r] = p->nz;
            int columns[9], count = 0;
            if (periodic && i == nx) {
                columns[count++] = r - nx;
                columns[count++] = r;
            } else
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++)
                        if (j + dy >= 0 && j + dy <= ny) {
                            int x = i + dx;
                            if (periodic) {
                                x = (x + nx) % nx;
                            } else if (x < 0 || x > nx)
                                continue;
                            columns[count++] = (j + dy) * (nx + 1) + x;
                        }
            /* ILU(0) requires strictly increasing columns in every CSR row. */
            for (int a = 1; a < count; a++) {
                int x = columns[a], b = a;
                while (b && columns[b - 1] > x) {
                    columns[b] = columns[b - 1];
                    b--;
                }
                columns[b] = x;
            }
            for (int a = 0; a < count; a++) {
                int c = columns[a];
                p->col[p->nz] = c;
                if (c == r)
                    p->diag[r] = p->nz;
                p->nz++;
            }
        }
    p->row[p->n] = p->nz;
    return 0;
}
/**
 * Create a reference transverse operator with nonperiodic longitude.
 */
int potential_create(Potential *p, int nx, int ny) {
    return create(p, nx, ny, 0);
}
/**
 * Create a reference transverse operator with explicit longitude periodicity.
 */
int potential_create_periodic(Potential *p, int nx, int ny) {
    return create(p, nx, ny, 1);
}
/**
 * Release sparse graph, factorization and iteration arrays.
 */
void potential_free(Potential *p) {
    free(p->row);
    free(p->col);
    free(p->diag);
    free(p->a);
    free(p->lu);
    free(p->rhs);
    free(p->phi);
    free(p->work);
    memset(p, 0, sizeof(*p));
}
/**
 * Clear current matrix/RHS values before assembling a new plasma state.
 */
void potential_zero(Potential *p) {
    memset(p->a, 0, p->nz * sizeof(double));
    memset(p->rhs, 0, p->n * sizeof(double));
}
/**
 * Locate and accumulate an entry in the preallocated sparse graph.
 */
static int entry(const Potential *p, int row, int col) {
    for (int k = p->row[row]; k < p->row[row + 1]; k++)
        if (p->col[k] == col)
            return k;
    return -1;
}
/**
 * Assemble one bilinear transverse element with both tensor cross terms and gravity-current
 * forcing.
 */
void potential_element(Potential *p, int i, int j, const double K[4], const double C[2]) {
    int id[4] = {j * (p->nx + 1) + i, j * (p->nx + 1) + i + 1, (j + 1) * (p->nx + 1) + i,
                 (j + 1) * (p->nx + 1) + i + 1};
    for (int a = 0; a < 4; a++)
        id[a] = canonical(p, id[a]);
    const double q[2] = {.5 - .5 / sqrt(3.), .5 + .5 / sqrt(3.)};
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++) {
            double grad[4][2];
            for (int k = 0; k < 4; k++) {
                int x = k % 2, y = k / 2;
                grad[k][0] = p->nx * (x ? 1 : -1) * (y ? q[b] : 1 - q[b]);
                grad[k][1] = p->ny * (y ? 1 : -1) * (x ? q[a] : 1 - q[a]);
            }
            for (int r = 0; r < 4; r++)
                if (id[r] >= p->first_row && id[r] < p->last_row && !boundary(p, id[r])) {
                    p->rhs[id[r]] += .25 * (grad[r][0] * C[0] + grad[r][1] * C[1]);
                    for (int c = 0; c < 4; c++)
                        if (!boundary(p, id[c]))
                            p->a[entry(p, id[r], id[c])] +=
                                .25 * (grad[r][0] * (K[0] * grad[c][0] + K[1] * grad[c][1]) +
                                       grad[r][1] * (K[2] * grad[c][0] + K[3] * grad[c][1]));
                }
        }
}
/**
 * Apply the reference solver's explicit physical/periodic potential boundary constraints.
 */
void potential_boundary(Potential *p) {
    for (int i = p->first_row; i < p->last_row; i++) {
        if (boundary(p, i)) {
            p->a[p->diag[i]] = 1;
            p->rhs[i] = p->phi[i] = 0;
        } else if (p->periodic_x && i % (p->nx + 1) == p->nx) {
            p->a[p->diag[i]] = 1;
            p->a[p->row[i]] = -1;
            p->rhs[i] = 0;
        }
    }
}
/**
 * Multiply the current sparse potential matrix by an input vector.
 */
static void apply(const Potential *p, const double *x, double *y) {
    for (int i = 0; i < p->n; i++) {
        y[i] = 0;
        for (int k = p->row[i]; k < p->row[i + 1]; k++)
            y[i] += p->a[k] * x[p->col[k]];
    }
}
/**
 * Evaluate a vector inner product for the reference iterative solver.
 */
static double dot(const double *a, const double *b, int n) {
    double s = 0;
    for (int i = 0; i < n; i++)
        s += a[i] * b[i];
    return s;
}
/**
 * Return the infinity norm of a vector.
 */
static double maxabs(const double *x, int n) {
    double m = 0;
    for (int i = 0; i < n; i++)
        m = fmax(m, fabs(x[i]));
    return m;
}
/**
 * Build the reference incomplete-LU preconditioner and reject unusable pivots.
 */
static int ilu(Potential *p) {
    memcpy(p->lu, p->a, p->nz * sizeof(double));
    for (int i = 0; i < p->n; i++) {
        for (int k = p->row[i]; k < p->diag[i]; k++) {
            int j = p->col[k];
            double pivot = p->lu[p->diag[j]];
            if (!isfinite(pivot) || pivot <= 0)
                return -1;
            double factor = p->lu[k] / pivot;
            p->lu[k] = factor;
            for (int z = p->diag[j] + 1; z < p->row[j + 1]; z++) {
                int target = entry(p, i, p->col[z]);
                if (target >= 0)
                    p->lu[target] -= factor * p->lu[z];
            }
        }
        if (!isfinite(p->lu[p->diag[i]]) || p->lu[p->diag[i]] <= 0)
            return -1;
    }
    return 0;
}
/**
 * Apply the stored incomplete-LU triangular factors.
 */
static void precondition(const Potential *p, const double *r, double *x) {
    for (int i = 0; i < p->n; i++) {
        x[i] = r[i];
        for (int k = p->row[i]; k < p->diag[i]; k++)
            x[i] -= p->lu[k] * x[p->col[k]];
    }
    for (int i = p->n - 1; i >= 0; i--) {
        for (int k = p->diag[i] + 1; k < p->row[i + 1]; k++)
            x[i] -= p->lu[k] * x[p->col[k]];
        x[i] /= p->lu[p->diag[i]];
    }
}
/**
 * Solve the assembled nonsymmetric system with preconditioned BiCGSTAB and a checked residual
 * criterion.
 */
int potential_solve(Potential *p, double rtol) {
    int n = p->n;
    p->iterations = 0;
    double scale = maxabs(p->rhs, n), tol = rtol * scale;
    if (scale == 0) {
        memset(p->phi, 0, n * sizeof(double));
        p->residual = 0;
        return 0;
    }
    if (ilu(p))
        return -1;
    double *r = p->work, *r0 = r + n, *v = r0 + n, *s = v + n, *t = s + n, *d = t + n, *dh = d + n,
           *sh = dh + n, *ax = sh + n;
    for (int restart = 0; restart < 3; restart++) {
        apply(p, p->phi, ax);
        for (int i = 0; i < n; i++)
            r[i] = p->rhs[i] - ax[i];
        if (maxabs(r, n) <= tol) {
            p->residual = maxabs(r, n) / scale;
            return p->iterations;
        }
        memcpy(r0, r, n * sizeof(double));
        memset(v, 0, n * sizeof(double));
        memset(d, 0, n * sizeof(double));
        double old = 1, alpha = 1, omega = 1;
        for (int iteration = 0; iteration < 2000; iteration++) {
            double rr = dot(r0, r, n);
            if (rr == 0 || !isfinite(rr) || omega == 0)
                break;
            double beta = (rr / old) * (alpha / omega);
            for (int i = 0; i < n; i++)
                d[i] = r[i] + beta * (d[i] - omega * v[i]);
            precondition(p, d, dh);
            apply(p, dh, v);
            double denom = dot(r0, v, n);
            if (denom == 0)
                break;
            alpha = rr / denom;
            for (int i = 0; i < n; i++)
                s[i] = r[i] - alpha * v[i];
            p->iterations++;
            if (maxabs(s, n) <= tol) {
                for (int i = 0; i < n; i++)
                    p->phi[i] += alpha * dh[i];
                break;
            }
            precondition(p, s, sh);
            apply(p, sh, t);
            double tt = dot(t, t, n);
            if (tt == 0)
                break;
            omega = dot(t, s, n) / tt;
            for (int i = 0; i < n; i++) {
                p->phi[i] += alpha * dh[i] + omega * sh[i];
                r[i] = s[i] - omega * t[i];
            }
            if (maxabs(r, n) <= tol)
                break;
            old = rr;
        }
    }
    apply(p, p->phi, ax);
    for (int i = 0; i < n; i++)
        r[i] = p->rhs[i] - ax[i];
    p->residual = maxabs(r, n) / scale;
    return p->residual <= 1.05 * rtol ? p->iterations : -2;
}
/**
 * Integrate local conductivity/current responses along field lines for the reference transverse
 * equation.
 */
void potential_coefficients(const TubeGrid *g, const double *rho, int first, int count, int hall,
                            double *a) {
    for (int p = first; p < first + count; p++) {
        double *v = a + 6 * p;
        memset(v, 0, 6 * sizeof(double));
        for (int k = 0; k < g->c.ns; k++) {
            int id = p * g->c.ns + k;
            const TubeCell *c = g->cell + id;
            double n = rho[id];
            double sp = c->pF * n + c->pE, sh = hall ? (c->hF * n + c->hE) : 0., w = c->volume;
            v[0] += w * sp * vdot(c->gx, c->gx);
            v[3] += w * sp * vdot(c->gy, c->gy);
            double symmetric = w * sp * vdot(c->gx, c->gy),
                   skew = w * sh * vdot(c->gx, vcross(c->b, c->gy));
            v[1] += symmetric + skew;
            v[2] += symmetric - skew;
            v[4] += w * (c->cF[0] * n + c->cE[0]);
            v[5] += w * (c->cF[1] * n + c->cE[1]);
        }
    }
}
/**
 * Assemble every transverse element and boundary condition from the current integrated
 * coefficients.
 */
void potential_assemble(Potential *p, const double *coefficients) {
    potential_zero(p);
    for (int j = 0; j < p->ny; j++)
        for (int i = 0; i < p->nx; i++)
            potential_element(p, i, j, coefficients + 6 * (j * p->nx + i),
                              coefficients + 6 * (j * p->nx + i) + 4);
    if (p->neumann_x) {
        /* Regional condition d_xi(phi_pol)=0, NOT zero total normal current.
         * From integration by parts, subtract the side K_xy*d_eta(phi)
         * bilinear term and C_x load. Extrapolate smooth buffer coefficients
         * to the side, retaining second-order boundary consistency. */
        for (int side = 0; side < 2; side++)
            for (int j = 0; j < p->ny; j++) {
                int i = side ? p->nx - 1 : 0, next = side ? i - 1 : i + 1, ix = side ? p->nx : 0,
                    normal = side ? 1 : -1;
                const double *edge = coefficients + 6 * (j * p->nx + i),
                             *inside = coefficients + 6 * (j * p->nx + next);
                double kxy = (1.5 * edge[1] - .5 * inside[1]) * p->nx * p->ny;
                double cx = (1.5 * edge[4] - .5 * inside[4]) * p->nx * p->ny;
                int id[2] = {j * (p->nx + 1) + ix, (j + 1) * (p->nx + 1) + ix};
                for (int r = 0; r < 2; r++)
                    if (id[r] >= p->first_row && id[r] < p->last_row && !boundary(p, id[r])) {
                        p->rhs[id[r]] -= normal * cx / (2 * p->ny);
                        for (int c = 0; c < 2; c++)
                            if (!boundary(p, id[c]))
                                p->a[entry(p, id[r], id[c])] -= normal * kxy * (c ? 1 : -1) * .5;
                    }
            }
    }
    potential_boundary(p);
}
