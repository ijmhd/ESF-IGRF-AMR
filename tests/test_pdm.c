/* SPDX-License-Identifier: GPL-3.0-only */
#include "tube.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "PDM geometry test failed line %d\n", __LINE__);                       \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)

/**
 * Return the mapped physical position used to build the synthetic reconstruction geometry.
 */
static V3 position(double x, double y, double z, int curved) {
    if (!curved)
        return (V3){x, y, z};
    return (V3){x * (1 + .25 * y + .15 * z), y * (1 + .2 * z), z};
}
/**
 * Evaluate its analytic Jacobian for independent quadrature/reconstruction checks.
 */
static double jacobian(double y, double z, int curved) {
    return curved ? (1 + .25 * y + .15 * z) * (1 + .2 * z) : 1;
}
/**
 * Evaluate the known polynomial used to test exact moment reproduction.
 */
static double polynomial(double x) {
    return 2 + .3 * x - .2 * x * x + .15 * x * x * x;
}
/**
 * Flatten fixture node indices with the field index contiguous.
 */
static int nid(const TubeGrid *g, int i, int j, int k) {
    return (j * (g->c.nx + 1) + i) * (g->c.ns + 1) + k;
}

/**
 * Build and test one mapped hexahedral geometry against polynomial and moment identities.
 */
static void geometry_test(int curved) {
    TubeGrid g = {.c = {.nx = 8, .ny = 8, .ns = 8}};
    int n = 8, sz = n * n * n;
    g.node = calloc((n + 1) * (n + 1) * (n + 1), sizeof(V3));
    g.cell = calloc(sz, sizeof(TubeCell));
    double *weights[3], *average[3];
    const double x[5] = {-.906179845938664, -.538469310105683, 0, .538469310105683,
                         .906179845938664};
    const double w[5] = {.236926885056189, .478628670499366, .568888888888889, .478628670499366,
                         .236926885056189};
    for (int j = 0; j <= n; j++)
        for (int i = 0; i <= n; i++)
            for (int k = 0; k <= n; k++)
                g.node[nid(&g, i, j, k)] =
                    position((double)i / n, (double)j / n, (double)k / n, curved);
    for (int a = 0; a < 3; a++) {
        weights[a] = calloc(4 * tube_face_count(&g, a), sizeof(double));
        average[a] = calloc(sz, sizeof(double));
    }
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++) {
                int id = tube_cell_id(&g, i, j, k);
                for (int a = 0; a < 5; a++)
                    for (int b = 0; b < 5; b++)
                        for (int c = 0; c < 5; c++) {
                            double q[3] = {(i + .5 + .5 * x[a]) / n, (j + .5 + .5 * x[b]) / n,
                                           (k + .5 + .5 * x[c]) / n};
                            double vol =
                                w[a] * w[b] * w[c] * jacobian(q[1], q[2], curved) / (8 * n * n * n);
                            g.cell[id].volume += vol;
                            for (int axis = 0; axis < 3; axis++)
                                average[axis][id] += vol * polynomial(q[axis]);
                        }
                for (int axis = 0; axis < 3; axis++)
                    average[axis][id] /= g.cell[id].volume;
            }
    CHECK(!tube_pdm_weights(&g, 0, weights));
    double error = 0, cartesian = 0;
    for (int axis = 0; axis < 3; axis++)
        for (int j = 0; j < n + (axis == 1); j++)
            for (int i = 0; i < n + (axis == 0); i++)
                for (int k = 0; k < n + (axis == 2); k++) {
                    int p[3] = {i, j, k}, fid = tube_face_id(&g, axis, i, j, k), f = p[axis];
                    if (f < 2 || f > n - 2) {
                        CHECK(isnan(weights[axis][4 * fid]));
                        continue;
                    }
                    double value = 0;
                    for (int c = 0; c < 4; c++) {
                        int cell[3] = {i, j, k};
                        cell[axis] += c - 2;
                        value += weights[axis][4 * fid + c] *
                                 average[axis][tube_cell_id(&g, cell[0], cell[1], cell[2])];
                        cartesian = fmax(cartesian,
                                         fabs(weights[axis][4 * fid + c] -
                                              (double[]){-1. / 12, 7. / 12, 7. / 12, -1. / 12}[c]));
                    }
                    error = fmax(error, fabs(value - polynomial((double)f / n)));
                }
    CHECK(error < 1e-12);
    if (!curved)
        CHECK(cartesian < 1e-12);
    printf(
        "PASS %s cubic physical-volume averages, max face error %.3g, Cartesian weight difference %.3g\n",
        curved ? "curved" : "Cartesian", error, cartesian);
    for (int a = 0; a < 3; a++) {
        free(weights[a]);
        free(average[a]);
    }
    tube_free(&g);
}

/**
 * Check mapped-cell quadrature, polynomial face reconstruction and limiter behavior.
 */
int main(void) {
    geometry_test(0);
    geometry_test(1);
    CHECK(tube_pdm_limit(2.5, 1, 2, 3, 1) == 2.5);
    CHECK(tube_pdm_limit(3, 1, 2, 3, .25) == 2.5);
    CHECK(tube_pdm_limit(3, 1, 2, 3, 0) == 2);
    CHECK(tube_pdm_limit(9, 1, 3, 2, 1) == 3);
    puts("PDM4 MOMENT AND LIMITER TESTS PASSED");
    return 0;
}
