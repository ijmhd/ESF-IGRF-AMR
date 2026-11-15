/*
 * Field-line-integrated electrostatic operator on constrained bilinear transverse elements. Both
 * Hall cross terms are retained; phi is constant along the traced field coordinate.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Q1 Galerkin assembly on canonical geometric subcells, constrained to the
 * active two-level conforming space. K contains the full Pedersen/Hall tensor.
 */
#include "amr3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.14159265358979323846
#define PCALL(call) PetscCallAbort(MPI_COMM_WORLD, (call))
typedef struct {
    int n;
    PetscInt id[16];
    PetscCount mo, ro;
    double *basis;
} CooElement;
typedef struct {
    int ne;
    PetscCount nm, nr;
    CooElement *e;
    PetscScalar *matrix, *rhs;
} Coo;
static void coo_create(A3 *a);
/**
 * Report whether a prescribed constant or time-dependent external eastward field is active.
 */
int a3_driver_active(const A3 *a) {
    return a->c.driver || a->c.background_east != 0;
}
/**
 * Evaluate the external eastward field at the actual SSPRK stage time.
 */
double a3_driver_field(const A3 *a) {
    double E = a->c.driver ? electric_driver_value(a->c.driver, a->time) : a->c.background_east;
    if (!isfinite(E))
        a3_fail((A3 *)a, "electric driver outside declared time interval");
    return E;
}
/**
 * Convert eastward field to the model's transverse potential-label gradient using the recorded
 * reference length.
 */
double a3_driver_gradient(const A3 *a) {
    return -a->east_scale_m * a3_driver_field(a);
}

/**
 * Flatten a canonical transverse potential vertex.
 */
static int vertex(A3 *a, int i, int j) {
    if (a->c.periodic_x && i == a->fx)
        i = 0;
    return j * (a->fx + 1) + i;
}
/**
 * Classify a potential node or apply its prescribed boundary role; inspect the explicit
 * regional/periodic branches below.
 */
static int boundary(A3 *a, int i, int j) {
    return j == 0 || j == a->fy || (!a->c.periodic_x && !a->c.regional_x && (i == 0 || i == a->fx));
}
/**
 * Accumulate one constrained-node basis contribution, merging duplicate independent DOFs.
 */
static void map_add(A3Map *map, int id, double w) {
    for (int k = 0; k < map->n; k++)
        if (map->id[k] == id) {
            map->w[k] += w;
            return;
        }
    if (map->n == 4)
        abort();
    map->id[map->n] = id;
    map->w[map->n++] = w;
}
/**
 * Create constrained field maps, distributed sparse PETSc objects, FGMRES/GAMG configuration and
 * optional precomputed COO basis tables.
 */
void a3_potential_create(A3 *a) {
    int nv = (a->fx + 1) * (a->fy + 1);
    a->map = a3_alloc(a, nv, sizeof(A3Map));
    a->phi = a3_alloc(a, nv, sizeof(double));
    if (a->c.amr_levels == 3)
        a3_macro_maps(a);
    else {
        int *slave = a3_alloc(a, nv, sizeof(int)), *dof = a3_alloc(a, nv, sizeof(int));
        for (int p = 0; p < nv; p++) {
            slave[p] = -1;
            dof[p] = -1;
        }
        for (int p = 0; p < a->nleaf; p++)
            if (!a->leaf[p].level) {
                A3Leaf l = a->leaf[p];
                for (int v = 0; v <= 2; v++)
                    for (int u = 0; u <= 2; u++)
                        if ((u & 1) || (v & 1))
                            slave[vertex(a, 2 * l.i + u, 2 * l.j + v)] = l.root;
            }
        for (int j = 0; j <= a->fy; j++)
            for (int i = 0; i < a->fx + !a->c.periodic_x; i++) {
                int id = vertex(a, i, j);
                if (!boundary(a, i, j) && slave[id] < 0)
                    dof[id] = a->ndof++;
            }
        for (int j = 0; j <= a->fy; j++)
            for (int i = 0; i < a->fx + !a->c.periodic_x; i++) {
                int id = vertex(a, i, j);
                if (boundary(a, i, j)) {
                    a->map[id].fixed = a->c.amr_case == 2 ? (double)j / a->fy : 0;
                    a->map[id].drive =
                        a3_driver_active(a) ? -a->east_scale_m * ((double)i / a->fx - .5) : 0;
                    continue;
                }
                if (slave[id] < 0)
                    map_add(a->map + id, dof[id], 1.);
            }
        for (int j = 0; j <= a->fy; j++)
            for (int i = 0; i < a->fx + !a->c.periodic_x; i++) {
                int id = vertex(a, i, j), b = slave[id];
                if (boundary(a, i, j) || b < 0)
                    continue;
                int x = 2 * (b % a->nx), y = 2 * (b / a->nx), fine_i = i;
                if (a->c.periodic_x && i == 0 && x == a->fx - 2)
                    fine_i = a->fx;
                double u = .5 * (fine_i - x), v = .5 * (j - y);
                for (int c = 0; c < 4; c++) {
                    double w = (c & 1 ? u : 1 - u) * (c >> 1 ? v : 1 - v);
                    A3Map *source = a->map + vertex(a, x + 2 * (c & 1), y + 2 * (c >> 1));
                    a->map[id].fixed += w * source->fixed;
                    a->map[id].drive += w * source->drive;
                    for (int k = 0; k < source->n; k++)
                        map_add(a->map + id, source->id[k], w * source->w[k]);
                }
            }
        if (a->c.periodic_x)
            for (int j = 0; j <= a->fy; j++)
                a->map[j * (a->fx + 1) + a->fx] = a->map[j * (a->fx + 1)];
        free(slave);
        free(dof);
    }
    int local_min = a->ndof / a->size, off_min = a->ndof - (a->ndof + a->size - 1) / a->size;
    PCALL(MatCreateAIJ(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, a->ndof, a->ndof,
                       local_min < 64 ? local_min : 64, NULL, off_min < 64 ? off_min : 64, NULL,
                       &a->matrix));
    PCALL(MatSetOption(a->matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    PCALL(MatCreateVecs(a->matrix, &a->solution, &a->rhs));
    PCALL(VecDuplicate(a->rhs, &a->residual));
    PCALL(VecSet(a->solution, 0));
    PCALL(VecScatterCreateToAll(a->solution, &a->scatter, &a->all));
    PCALL(KSPCreate(MPI_COMM_WORLD, &a->ksp));
    PCALL(KSPSetOperators(a->ksp, a->matrix, a->matrix));
    PCALL(KSPSetType(a->ksp, KSPFGMRES));
    PCALL(KSPGMRESSetRestart(a->ksp, 40));
    PCALL(KSPSetPCSide(a->ksp, PC_RIGHT));
    PCALL(KSPSetInitialGuessNonzero(a->ksp, PETSC_TRUE));
    PC pc;
    PCALL(KSPGetPC(a->ksp, &pc));
    PCALL(PCSetType(pc, PCGAMG));
    PCALL(PCGAMGSetNSmooths(pc, 1));
    PCALL(PCGAMGSetReuseInterpolation(pc, PETSC_TRUE));
    PCALL(KSPSetNormType(a->ksp, KSP_NORM_UNPRECONDITIONED));
    if (a->c.amr_coo)
        coo_create(a);
}
/**
 * Precompute constrained bilinear element bases and PETSc COO sparsity/index maps; material
 * coefficients remain current-stage data.
 */
static void coo_create(A3 *a) {
    Coo *coo = a3_alloc(a, 1, sizeof(Coo));
    a->coo_context = coo;
    coo->ne = 2 * (a->row_last - a->row_first) * a->fx;
    coo->e = a3_alloc(a, coo->ne, sizeof(CooElement));
    double B[4][4][4] = {0}, L[2][4] = {0};
    const double q[2] = {.5 - .5 / sqrt(3.), .5 + .5 / sqrt(3.)};
    for (int v = 0; v < 2; v++)
        for (int u = 0; u < 2; u++) {
            double g[4][2];
            for (int c = 0; c < 4; c++) {
                g[c][0] = a->fx * (c & 1 ? 1 : -1) * (c >> 1 ? q[v] : 1 - q[v]);
                g[c][1] = a->fy * (c >> 1 ? 1 : -1) * (c & 1 ? q[u] : 1 - q[u]);
            }
            for (int r = 0; r < 4; r++) {
                for (int d = 0; d < 2; d++)
                    L[d][r] += .25 * g[r][d];
                for (int c = 0; c < 4; c++)
                    for (int d = 0; d < 4; d++)
                        B[d][r][c] += .25 * g[r][d / 2] * g[c][d % 2];
            }
        }
    int elem = 0;
    for (int j = 2 * a->row_first; j < 2 * a->row_last; j++)
        for (int i = 0; i < a->fx; i++) {
            CooElement *e = coo->e + elem++;
            A3Map *map[4] = {a->map + vertex(a, i, j), a->map + vertex(a, i + 1, j),
                             a->map + vertex(a, i, j + 1), a->map + vertex(a, i + 1, j + 1)};
            int index[4][4];
            for (int c = 0; c < 4; c++)
                for (int k = 0; k < map[c]->n; k++) {
                    int position = 0;
                    while (position < e->n && e->id[position] != map[c]->id[k])
                        position++;
                    if (position == e->n)
                        e->id[e->n++] = map[c]->id[k];
                    index[c][k] = position;
                }
            e->mo = coo->nm;
            e->ro = coo->nr;
            coo->nm += (PetscCount)e->n * e->n;
            coo->nr += e->n;
            int n = e->n, nn = n * n;
            e->basis = a3_alloc(a, (size_t)4 * nn + 10 * n, sizeof(double));
            double *rhs = e->basis + 4 * nn;
            for (int r = 0; r < 4; r++)
                for (int rr = 0; rr < map[r]->n; rr++) {
                    int R = index[r][rr];
                    double w = map[r]->w[rr];
                    for (int d = 0; d < 2; d++)
                        rhs[d * n + R] += w * L[d][r];
                    for (int c = 0; c < 4; c++)
                        for (int d = 0; d < 4; d++) {
                            rhs[(2 + d) * n + R] -= w * B[d][r][c] * map[c]->fixed;
                            rhs[(6 + d) * n + R] -= w * B[d][r][c] * map[c]->drive;
                            for (int cc = 0; cc < map[c]->n; cc++)
                                e->basis[d * nn + R * n + index[c][cc]] +=
                                    w * B[d][r][c] * map[c]->w[cc];
                        }
                }
        }
    PetscCount matrix_count = coo->nm, rhs_count = coo->nr;
    if (a->c.regional_x && !a->c.periodic_x)
        for (int side = 0; side < 2; side++)
            for (int j = 2 * a->row_first; j < 2 * a->row_last; j++) {
                int x = side ? a->fx : 0,
                    n = a->map[vertex(a, x, j)].n + a->map[vertex(a, x, j + 1)].n;
                matrix_count += (PetscCount)n * n;
                rhs_count += n;
            }
    PetscInt *ri = a3_alloc(a, (size_t)matrix_count, sizeof(PetscInt)),
             *ci = a3_alloc(a, (size_t)matrix_count, sizeof(PetscInt)),
             *vi = a3_alloc(a, (size_t)rhs_count, sizeof(PetscInt));
    for (int p = 0; p < coo->ne; p++) {
        CooElement *e = coo->e + p;
        for (int r = 0; r < e->n; r++) {
            vi[e->ro + r] = e->id[r];
            for (int c = 0; c < e->n; c++) {
                PetscCount k = e->mo + r * e->n + c;
                ri[k] = e->id[r];
                ci[k] = e->id[c];
            }
        }
    }
    PetscCount mi = coo->nm, vr = coo->nr;
    if (a->c.regional_x && !a->c.periodic_x)
        for (int side = 0; side < 2; side++)
            for (int j = 2 * a->row_first; j < 2 * a->row_last; j++) {
                int x = side ? a->fx : 0;
                A3Map *map[2] = {a->map + vertex(a, x, j), a->map + vertex(a, x, j + 1)};
                for (int r = 0; r < 2; r++)
                    for (int rr = 0; rr < map[r]->n; rr++) {
                        vi[vr++] = map[r]->id[rr];
                        for (int c = 0; c < 2; c++)
                            for (int cc = 0; cc < map[c]->n; cc++) {
                                ri[mi] = map[r]->id[rr];
                                ci[mi++] = map[c]->id[cc];
                            }
                    }
            }
    if (mi != matrix_count || vr != rhs_count)
        a3_fail(a, "COO layout counts");
    PCALL(MatSetPreallocationCOO(a->matrix, matrix_count, ri, ci));
    PCALL(VecSetPreallocationCOO(a->rhs, rhs_count, vi));
    free(ri);
    free(ci);
    free(vi);
    coo->matrix = a3_alloc(a, (size_t)matrix_count, sizeof(PetscScalar));
    coo->rhs = a3_alloc(a, (size_t)rhs_count, sizeof(PetscScalar));
}
/**
 * Integrate both Pedersen/Hall tensor cross terms and gravity-current forcing along the full
 * field line using canonical physical volumes.
 */
static void integrate_tensor(A3 *a, int i, int j, double K[4], double C[2]) {
    memset(K, 0, 4 * sizeof(double));
    memset(C, 0, 2 * sizeof(double));
    if (a->c.amr_case == 2) {
        double x = (i + .5) / a->fx, y = (j + .5) / a->fy, area = 1. / (a->fx * a->fy);
        K[0] = (2 + .2 * cos(2 * PI * x)) * area;
        K[3] = (1.5 + .1 * y) * area;
        K[1] = (.25 + .12 * sin(2 * PI * x)) * area;
        K[2] = (.25 - .12 * sin(2 * PI * x)) * area;
        C[0] = K[1];
        C[1] = K[3];
        double p = a3_driver_gradient(a);
        C[0] += K[0] * p;
        C[1] += K[2] * p;
        return;
    }
    for (int k = 0; k < a->ns; k++) {
        int id = tube_cell_id(a->g, i, j, k);
        double *r = a->response + 12 * id, n = a3_sample(a, 0, 1, i, j, k);
        for (int d = 0; d < 4; d++)
            K[d] += r[d] * n + r[4 + d];
        C[0] += r[8] * n + r[10];
        C[1] += r[9] * n + r[11];
    }
}
/**
 * Seed independent potential DOFs from a prior canonical potential for a warm start after
 * regridding.
 */
void a3_potential_seed(A3 *a, const double *phi) {
    PetscInt first, last;
    PCALL(VecGetOwnershipRange(a->solution, &first, &last));
    for (int v = 0; v < (a->fx + 1) * (a->fy + 1); v++) {
        A3Map *p = a->map + v;
        if (p->n == 1 && p->w[0] == 1 && p->id[0] >= first && p->id[0] < last)
            PCALL(VecSetValue(a->solution, p->id[0], phi[v], INSERT_VALUES));
    }
    PCALL(VecAssemblyBegin(a->solution));
    PCALL(VecAssemblyEnd(a->solution));
}
/**
 * Assemble one bilinear element with the full nonsymmetric tensor and prescribed potential
 * contributions using the reference element path.
 */
static void element(A3 *a, int i, int j, const double K[4], const double C[2]) {
    int id[4] = {vertex(a, i, j), vertex(a, i + 1, j), vertex(a, i, j + 1),
                 vertex(a, i + 1, j + 1)};
    const double q[2] = {.5 - .5 / sqrt(3.), .5 + .5 / sqrt(3.)};
    double matrix[4][4] = {{0}}, load[4] = {0};
    for (int v = 0; v < 2; v++)
        for (int u = 0; u < 2; u++) {
            double grad[4][2];
            for (int c = 0; c < 4; c++) {
                grad[c][0] = a->fx * (c & 1 ? 1 : -1) * (c >> 1 ? q[v] : 1 - q[v]);
                grad[c][1] = a->fy * (c >> 1 ? 1 : -1) * (c & 1 ? q[u] : 1 - q[u]);
            }
            for (int r = 0; r < 4; r++) {
                load[r] += .25 * (grad[r][0] * C[0] + grad[r][1] * C[1]);
                for (int c = 0; c < 4; c++)
                    matrix[r][c] += .25 * (grad[r][0] * (K[0] * grad[c][0] + K[1] * grad[c][1]) +
                                           grad[r][1] * (K[2] * grad[c][0] + K[3] * grad[c][1]));
            }
        }
    PetscInt ids[16];
    int count = 0, index[4][4];
    for (int c = 0; c < 4; c++) {
        A3Map *map = a->map + id[c];
        for (int k = 0; k < map->n; k++) {
            int position = 0;
            while (position < count && ids[position] != map->id[k])
                position++;
            if (position == count)
                ids[count++] = map->id[k];
            index[c][k] = position;
        }
    }
    if (!count)
        return;
    double block[256] = {0}, rhs[16] = {0};
    for (int r = 0; r < 4; r++) {
        A3Map *row = a->map + id[r];
        double value = load[r];
        for (int c = 0; c < 4; c++)
            value -=
                matrix[r][c] * (a->map[id[c]].fixed + a3_driver_field(a) * a->map[id[c]].drive);
        for (int rr = 0; rr < row->n; rr++) {
            int R = index[r][rr];
            rhs[R] += row->w[rr] * value;
            for (int c = 0; c < 4; c++) {
                A3Map *col = a->map + id[c];
                for (int cc = 0; cc < col->n; cc++)
                    block[R * count + index[c][cc]] += row->w[rr] * matrix[r][c] * col->w[cc];
            }
        }
    }
    PCALL(VecSetValues(a->rhs, count, ids, rhs, ADD_VALUES));
    PCALL(MatSetValues(a->matrix, count, ids, count, ids, block, ADD_VALUES));
}
/**
 * Impose the regional polarization-potential normal-derivative relation, retaining Hall and
 * nonorthogonal metric terms.
 */
static void regional_sides(A3 *a) {
    if (a->c.periodic_x || !a->c.regional_x)
        return;
    /* d_xi(phi_total)=d_xi(phi_ext), i.e. d_xi(phi_pol)=0. Retain
     * tangential Hall/geometry current and constrain hanging boundary nodes. */
    for (int side = 0; side < 2; side++)
        for (int j = 2 * a->row_first; j < 2 * a->row_last; j++) {
            int i = side ? a->fx - 1 : 0, next = side ? i - 1 : i + 1, ix = side ? a->fx : 0,
                normal = side ? 1 : -1;
            double K[4], C[2], innerK[4], innerC[2];
            integrate_tensor(a, i, j, K, C);
            integrate_tensor(a, next, j, innerK, innerC);
            double kxy = (1.5 * K[1] - .5 * innerK[1]) * a->fx * a->fy;
            double cx = (1.5 * C[0] - .5 * innerC[0]) * a->fx * a->fy;
            double kxx = (1.5 * K[0] - .5 * innerK[0]) * a->fx * a->fy;
            A3Map *map[2] = {a->map + vertex(a, ix, j), a->map + vertex(a, ix, j + 1)};
            for (int r = 0; r < 2; r++)
                for (int rr = 0; rr < map[r]->n; rr++) {
                    PetscInt R = map[r]->id[rr];
                    double load = normal * (kxx * a3_driver_gradient(a) - cx) / (2 * a->fy);
                    for (int c = 0; c < 2; c++)
                        load -= (-normal * kxy * (c ? 1 : -1) * .5) *
                                (map[c]->fixed + a3_driver_field(a) * map[c]->drive);
                    load *= map[r]->w[rr];
                    PCALL(VecSetValue(a->rhs, R, load, ADD_VALUES));
                    for (int c = 0; c < 2; c++)
                        for (int cc = 0; cc < map[c]->n; cc++) {
                            PetscInt Col = map[c]->id[cc];
                            double value =
                                map[r]->w[rr] * (-normal * kxy * (c ? 1 : -1) * .5) * map[c]->w[cc];
                            PCALL(MatSetValue(a->matrix, R, Col, value, ADD_VALUES));
                        }
                }
        }
}
/**
 * Assemble the current-stage integrated operator, solve with controlled preconditioner reuse,
 * explicitly verify the true residual, and scatter canonical potential. Collective; no worker-
 * thread PETSc calls.
 */
void a3_potential(A3 *a) {
    double start = wall3();
    PCALL(MatZeroEntries(a->matrix));
    PCALL(VecSet(a->rhs, 0));
    Coo *coo = a->coo_context;
    int elem = 0;
    for (int j = 2 * a->row_first; j < 2 * a->row_last; j++)
        for (int i = 0; i < a->fx; i++) {
            double K[4], C[2];
            integrate_tensor(a, i, j, K, C);
            if (!coo)
                element(a, i, j, K, C);
            else {
                CooElement *e = coo->e + elem++;
                int n = e->n, nn = n * n;
                double *b = e->basis, *rhs = b + 4 * nn;
                for (int k = 0; k < nn; k++)
                    coo->matrix[e->mo + k] = b[k] * K[0] + b[nn + k] * K[1] + b[2 * nn + k] * K[2] +
                                             b[3 * nn + k] * K[3];
                for (int k = 0; k < n; k++)
                    coo->rhs[e->ro + k] =
                        rhs[k] * C[0] + rhs[n + k] * C[1] + rhs[2 * n + k] * K[0] +
                        rhs[3 * n + k] * K[1] + rhs[4 * n + k] * K[2] + rhs[5 * n + k] * K[3] +
                        a3_driver_field(a) * (rhs[6 * n + k] * K[0] + rhs[7 * n + k] * K[1] +
                                              rhs[8 * n + k] * K[2] + rhs[9 * n + k] * K[3]);
            }
        }
    a->t_coeff += wall3() - start;
    start = wall3();
    if (coo) {
        PCALL(MatSetValuesCOO(a->matrix, coo->matrix, INSERT_VALUES));
        PCALL(VecSetValuesCOO(a->rhs, coo->rhs, INSERT_VALUES));
    }
    regional_sides(a);
    PCALL(MatAssemblyBegin(a->matrix, MAT_FINAL_ASSEMBLY));
    PCALL(MatAssemblyEnd(a->matrix, MAT_FINAL_ASSEMBLY));
    PCALL(VecAssemblyBegin(a->rhs));
    PCALL(VecAssemblyEnd(a->rhs));
    PetscReal scale;
    PCALL(VecNorm(a->rhs, NORM_INFINITY, &scale));
    a->potential_residual = 0;
    int iterations = 0;
    if (scale == 0)
        PCALL(VecSet(a->solution, 0));
    else
        for (int attempt = 0; attempt < 2; attempt++) {
            PCALL(KSPSetReusePreconditioner(
                a->ksp,
                !attempt && a->solves % 32 && a->last_iterations < 20 ? PETSC_TRUE : PETSC_FALSE));
            PCALL(KSPSetTolerances(a->ksp, 1e-14, .25 * a->c.rtol * scale, PETSC_DEFAULT, 400));
            PCALL(KSPSolve(a->ksp, a->rhs, a->solution));
            PetscInt count;
            PCALL(KSPGetIterationNumber(a->ksp, &count));
            iterations += count;
            PCALL(MatMult(a->matrix, a->solution, a->residual));
            PCALL(VecAYPX(a->residual, -1, a->rhs));
            PetscReal norm;
            PCALL(VecNorm(a->residual, NORM_INFINITY, &norm));
            a->potential_residual = norm / scale;
            if (isfinite(a->potential_residual) && a->potential_residual <= 1.05 * a->c.rtol)
                break;
        }
    if (!isfinite(a->potential_residual) || a->potential_residual > 1.05 * a->c.rtol)
        a3_fail(a, "composite tensor potential failed true residual");
    a->last_iterations = iterations;
    a->solves++;
    a->max_residual = fmax(a->max_residual, a->potential_residual);
    PCALL(VecScatterBegin(a->scatter, a->solution, a->all, INSERT_VALUES, SCATTER_FORWARD));
    PCALL(VecScatterEnd(a->scatter, a->solution, a->all, INSERT_VALUES, SCATTER_FORWARD));
    const PetscScalar *phi;
    PCALL(VecGetArrayRead(a->all, &phi));
    for (int v = 0; v < (a->fx + 1) * (a->fy + 1); v++) {
        A3Map *map = a->map + v;
        a->phi[v] = map->fixed + a3_driver_field(a) * map->drive;
        for (int k = 0; k < map->n; k++)
            a->phi[v] += map->w[k] * phi[map->id[k]];
    }
    PCALL(VecRestoreArrayRead(a->all, &phi));
    a->t_potential += wall3() - start;
}
/**
 * Destroy constrained maps, COO workspaces and every PETSc object owned by this mesh.
 */
void a3_potential_free(A3 *a) {
    if (a->ksp) {
        PCALL(VecScatterDestroy(&a->scatter));
        PCALL(VecDestroy(&a->all));
        PCALL(VecDestroy(&a->rhs));
        PCALL(VecDestroy(&a->solution));
        PCALL(VecDestroy(&a->residual));
        PCALL(KSPDestroy(&a->ksp));
        PCALL(MatDestroy(&a->matrix));
    }
    Coo *coo = a->coo_context;
    if (coo) {
        for (int p = 0; p < coo->ne; p++)
            free(coo->e[p].basis);
        free(coo->e);
        free(coo->matrix);
        free(coo->rhs);
        free(coo);
        a->coo_context = NULL;
    }
    free(a->map);
    free(a->phi);
}
