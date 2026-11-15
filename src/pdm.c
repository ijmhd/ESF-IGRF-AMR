/*
 * Physical-volume four-cell cubic reconstruction and partial-donor-cell limiting. PDM4 describes
 * the spatial stencil; the coupled time update remains second order.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Physical-volume moments for a four-cell cubic face reconstruction.
 * PDM clipping follows the workspace's audited GPLv3 C kernel and PDMU7.m
 * in bzhang1982/gamera_matlab. See THIRD_PARTY.md.
 */
#include "tube.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * Flatten node coordinates with the along-field coordinate contiguous.
 */
static int node_id(const TubeGrid *g, int i, int j, int k) {
    return ((j * (g->c.nx + 1) + i) * (g->c.ns + 1) + k);
}

/**
 * Evaluate the signed Jacobian determinant of a trilinear hexahedron at one logical quadrature
 * point.
 */
static double determinant(const V3 v[8], const double q[3]) {
    V3 e[3] = {{0, 0, 0}};
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++)
            for (int c = 0; c < 2; c++) {
                V3 p = v[(a * 2 + b) * 2 + c];
                double x = a ? q[0] : 1 - q[0], y = b ? q[1] : 1 - q[1], z = c ? q[2] : 1 - q[2];
                e[0] = vadd(e[0], vmul(p, (a ? 1 : -1) * y * z));
                e[1] = vadd(e[1], vmul(p, x * (b ? 1 : -1) * z));
                e[2] = vadd(e[2], vmul(p, x * y * (c ? 1 : -1)));
            }
    return vdot(e[0], vcross(e[1], e[2]));
}

/**
 * Integrate normalized physical-volume moments up to cubic degree using tensor three-point
 * Gaussian quadrature; reject folded cells and inconsistent volumes.
 */
int tube_cell_moments(const TubeGrid *g, int i, int j, int k, double out[12]) {
    const double point[] = {.5 - sqrt(3. / 5) / 2, .5, .5 + sqrt(3. / 5) / 2};
    const double weight[] = {5. / 18, 4. / 9, 5. / 18};
    V3 v[8];
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++)
            for (int c = 0; c < 2; c++)
                v[(a * 2 + b) * 2 + c] = g->node[node_id(g, i + a, j + b, k + c)];
    memset(out, 0, 12 * sizeof(double));
    double sign = determinant(v, (double[]){.5, .5, .5});
    if (!isfinite(sign) || sign == 0)
        return -1;
    /* det(J) is at most quadratic in each trilinear label. Three-point
     * quadrature integrates det(J)*q^p exactly for p <= 3. */
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            for (int c = 0; c < 3; c++) {
                double q[] = {point[a], point[b], point[c]}, jac = determinant(v, q);
                if (!isfinite(jac) || jac * sign <= 0)
                    return -1;
                double w = fabs(jac) * weight[a] * weight[b] * weight[c];
                for (int axis = 0; axis < 3; axis++) {
                    double power = 1;
                    for (int degree = 0; degree < 4; degree++) {
                        out[4 * axis + degree] += w * power;
                        power *= q[axis];
                    }
                }
            }
    double vol = g->cell[tube_cell_id(g, i, j, k)].volume;
    if (!(vol > 0) || fabs(out[0] / vol - 1) > 1e-9)
        return -1;
    for (int axis = 0; axis < 3; axis++)
        for (int degree = 0; degree < 4; degree++)
            out[4 * axis + degree] /= vol;
    return 0;
}

/**
 * Solve the four-by-four moment matching system for a face value, using pivoted elimination and
 * polynomial reproduction checks.
 */
int tube_moment_weights(const double moments[4][4], double weights[4]) {
    double a[4][5];
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++)
            a[r][c] = moments[c][r];
        a[r][4] = r == 0;
    }
    for (int c = 0; c < 4; c++) {
        int pivot = c;
        for (int r = c + 1; r < 4; r++)
            if (fabs(a[r][c]) > fabs(a[pivot][c]))
                pivot = r;
        if (fabs(a[pivot][c]) < 1e-12)
            return -1;
        for (int k = c; k < 5; k++) {
            double tmp = a[c][k];
            a[c][k] = a[pivot][k];
            a[pivot][k] = tmp;
        }
        double diagonal = a[c][c];
        for (int k = c; k < 5; k++)
            a[c][k] /= diagonal;
        for (int r = 0; r < 4; r++)
            if (r != c) {
                double factor = a[r][c];
                for (int k = c; k < 5; k++)
                    a[r][k] -= factor * a[c][k];
            }
    }
    for (int c = 0; c < 4; c++) {
        weights[c] = a[c][4];
        if (!isfinite(weights[c]))
            return -1;
    }
    for (int degree = 0; degree < 4; degree++) {
        double value = 0;
        for (int c = 0; c < 4; c++)
            value += weights[c] * moments[c][degree];
        if (fabs(value - (degree == 0)) > 1e-10)
            return -1;
    }
    return 0;
}

/**
 * Prepare four-cell cubic face weights from physical moments for a uniform logical grid;
 * unavailable boundary stencils remain NaN and use a lower-order branch.
 */
int tube_pdm_weights(const TubeGrid *g, int periodic_x, double *weights[3]) {
    int nx = g->c.nx, ny = g->c.ny, ns = g->c.ns, n = nx * ny * ns;
    if (!g->node)
        return -1;
    double *moment = calloc((size_t)n * 12, sizeof(double));
    if (!moment)
        return -1;
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            for (int k = 0; k < ns; k++) {
                if (tube_cell_moments(g, i, j, k, moment + 12 * tube_cell_id(g, i, j, k))) {
                    free(moment);
                    return -1;
                }
            }
    const int binomial[4][4] = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 2, 1, 0}, {1, 3, 3, 1}};
    int count[3] = {nx, ny, ns};
    for (int axis = 0; axis < 3; axis++) {
        for (int f = 0; f < 4 * tube_face_count(g, axis); f++)
            weights[axis][f] = NAN;
        for (int j = 0; j < ny + (axis == 1); j++)
            for (int i = 0; i < nx + (axis == 0); i++)
                for (int k = 0; k < ns + (axis == 2); k++) {
                    int p[] = {i, j, k}, face = p[axis];
                    if (!(axis == 0 && periodic_x) && (face < 2 || face > count[axis] - 2))
                        continue;
                    double moments[4][4] = {{0}};
                    for (int c = 0; c < 4; c++) {
                        int cell[3] = {i, j, k}, offset = c - 2;
                        cell[axis] = face + offset;
                        if (axis == 0 && periodic_x)
                            cell[0] = (cell[0] + nx) % nx;
                        const double *local =
                            moment + 12 * tube_cell_id(g, cell[0], cell[1], cell[2]) + 4 * axis;
                        for (int degree = 0; degree < 4; degree++)
                            for (int r = 0; r <= degree; r++)
                                moments[c][degree] += binomial[degree][r] *
                                                      pow((double)offset, degree - r) * local[r];
                    }
                    if (tube_moment_weights((const double(*)[4])moments,
                                            weights[axis] + 4 * tube_face_id(g, axis, i, j, k))) {
                        free(moment);
                        return -1;
                    }
                }
    }
    free(moment);
    return 0;
}

/**
 * Clamp the high-order state to donor/downstream bounds, suppress extrema and cap the correction
 * by 2*beta times the upstream difference.
 */
double tube_pdm_limit(double high, double upstream, double donor, double downstream, double beta) {
    double a = donor - upstream, b = downstream - donor;
    double clipped = fmax(fmin(donor, downstream), fmin(fmax(donor, downstream), high));
    if (a * b <= 0)
        return donor;
    return donor + copysign(fmin(fabs(clipped - donor), 2 * beta * fabs(a)), b);
}
