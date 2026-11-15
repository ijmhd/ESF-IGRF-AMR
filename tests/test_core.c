#include "potential.h"
#include "integrator.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PI 3.14159265358979323846
#define CHECK(c)                                                                                   \
    do {                                                                                           \
        if (!(c)) {                                                                                \
            fprintf(stderr, "check failed %s:%d: %s\n", __FILE__, __LINE__, #c);                   \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
/**
 * Run the basic manufactured-potential verification and report its error.
 */
static double mms(int n) {
    Potential p;
    CHECK(!potential_create(&p, n, n));
    potential_zero(&p);
    const double q[2] = {.5 - .5 / sqrt(3.), .5 + .5 / sqrt(3.)};
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            double x = (i + .5) / n, y = (j + .5) / n,
                   K[4] = {(2 + .2 * x) / (n * n), (.35 + .1 * y) / (n * n),
                           (.05 - .1 * y) / (n * n), (1.5 + .3 * y) / (n * n)},
                   C[2] = {0, 0};
            potential_element(&p, i, j, K, C);
            for (int a = 0; a < 2; a++)
                for (int b = 0; b < 2; b++) {
                    x = (i + q[a]) / n;
                    y = (j + q[b]) / n;
                    double f = PI * PI * (3.5 + .2 * x + .3 * y) * sin(PI * x) * sin(PI * y) -
                               .4 * PI * PI * cos(PI * x) * cos(PI * y) -
                               .1 * PI * cos(PI * x) * sin(PI * y) -
                               .3 * PI * sin(PI * x) * cos(PI * y);
                    for (int k = 0; k < 4; k++) {
                        int ix = i + k % 2, iy = j + k / 2;
                        if (ix && ix < n && iy && iy < n)
                            p.rhs[iy * (n + 1) + ix] += .25 / (n * n) * f *
                                                        (k % 2 ? q[a] : 1 - q[a]) *
                                                        (k / 2 ? q[b] : 1 - q[b]);
                    }
                }
        }
    potential_boundary(&p);
    CHECK(potential_solve(&p, 1e-11) >= 0);
    double e = 0;
    for (int j = 1; j < n; j++)
        for (int i = 1; i < n; i++)
            e += pow(p.phi[j * (n + 1) + i] - sin(PI * i / n) * sin(PI * j / n), 2);
    e = sqrt(e / (n * n));
    printf("tensor/Hall MMS n=%d error=%.12g residual=%.3g\n", n, e, p.residual);
    potential_free(&p);
    return e;
}
/**
 * Integrate the test electric-field/flux diagnostic for a manufactured geometry.
 */
static double E_integral(const Igrf *m, int n) {
    V3 *line = calloc(n + 1, sizeof(V3));
    CHECK(!tube_trace(m, 114.69267876, 900e3, 90e3, 2000, n, line));
    TubeGrid g = {0};
    g.c.nx = g.c.ny = 1;
    g.c.ns = n;
    g.cell = calloc(n, sizeof(TubeCell));
    for (int k = 0; k < n; k++) {
        TubeCell *c = g.cell + k;
        c->position = vmul(vadd(line[k], line[k + 1]), .5);
        c->height = vnorm(c->position) - 6371200.;
        V3 b = igrf_field(m, c->position);
        c->B = vnorm(b);
        c->b = vunit(b);
    }
    tube_conductivity(&g, 1e10);
    double sigma = 0;
    for (int k = 0; k < n; k++)
        sigma += g.cell[k].pE * vnorm(vsub(line[k + 1], line[k]));
    free(g.cell);
    free(line);
    return sigma;
}
/**
 * Run the regional-side manufactured solution, retaining the relevant cross terms.
 */
static double regional_mms(int n) {
    Potential p;
    CHECK(!potential_create(&p, n, n));
    p.neumann_x = 1;
    double *coeff = calloc(6 * n * n, sizeof(double));
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            double x = (i + .5) / n, y = (j + .5) / n, *a = coeff + 6 * (j * n + i),
                   area = 1. / (n * n);
            double dx = -PI * sin(PI * x) * sin(PI * y), dy = PI * cos(PI * x) * cos(PI * y);
            a[0] = (2 + .2 * x) * area;
            a[1] = (.35 + .1 * y) * area;
            a[2] = (.05 - .1 * y) * area;
            a[3] = (1.5 + .3 * y) * area;
            /* Analytic current J=(cos(pi*y),0) crosses the side boundaries.
             * Zero-normal-current boundaries would not solve this problem. */
            a[4] = a[0] * dx + a[1] * dy + area * cos(PI * y);
            a[5] = a[2] * dx + a[3] * dy;
        }
    potential_assemble(&p, coeff);
    CHECK(potential_solve(&p, 1e-11) >= 0);
    double error = 0;
    for (int j = 1; j < n; j++)
        for (int i = 0; i <= n; i++)
            error += pow(p.phi[j * (n + 1) + i] - cos(PI * i / n) * sin(PI * j / n), 2);
    error = sqrt(error / (n * n));
    printf("regional nonorthogonal/Hall boundary MMS %d error %.12g residual %.3g\n", n, error,
           p.residual);
    free(coeff);
    potential_free(&p);
    return error;
}
/**
 * Run the periodic-longitude manufactured solution and check its independent-node constraints.
 */
static double periodic_mms(int n) {
    Potential p;
    CHECK(!potential_create_periodic(&p, n, n));
    double *coeff = calloc(6 * n * n, sizeof(double));
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            double x = (i + .5) / n, y = (j + .5) / n, *a = coeff + 6 * (j * n + i),
                   area = 1. / (n * n);
            double dx = 2 * PI * cos(2 * PI * x) * sin(PI * y),
                   dy = PI * sin(2 * PI * x) * cos(PI * y);
            a[0] = (2 + .2 * sin(2 * PI * x)) * area;
            a[1] = (.35 + .1 * y) * area;
            a[2] = (.05 - .1 * y) * area;
            a[3] = (1.5 + .3 * y) * area;
            a[4] = a[0] * dx + a[1] * dy;
            a[5] = a[2] * dx + a[3] * dy;
        }
    potential_assemble(&p, coeff);
    CHECK(potential_solve(&p, 1e-11) >= 0);
    double error = 0;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < n; i++)
            error += pow(p.phi[j * (n + 1) + i] - sin(2 * PI * i / n) * sin(PI * j / n), 2);
    for (int j = 0; j <= n; j++)
        CHECK(fabs(p.phi[j * (n + 1)] - p.phi[j * (n + 1) + n]) < 1e-8);
    error = sqrt(error / (n * n));
    printf("periodic nonorthogonal/Hall MMS %d error %.12g residual %.3g\n", n, error, p.residual);
    free(coeff);
    potential_free(&p);
    return error;
}
/**
 * Check magnetic geometry and manufactured transverse potential problems, including regional and
 * periodic boundaries.
 */
int main(int argc, char **argv) {
    for (int i = -12; i <= 12; i++) {
        double nu = pow(10., i), decay, p1, p2;
        exponential_weights(nu, 2., &decay, &p1, &p2);
        CHECK(p1 >= p2 && p2 >= 0);
        double updated = decay * .3 + 2 * ((p1 - p2) * .7 + p2 * .7);
        double exact = decay * .3 + .7 * (-expm1(-nu * 2)) / nu;
        CHECK(fabs(updated - exact) < 2e-14 * fmax(1., fabs(exact)));
    }
    CHECK(argc == 2);
    Igrf m;
    CHECK(!igrf_load(&m, argv[1], 2025, 1));
    V3 line[65];
    CHECK(!tube_trace(&m, 114.2, 900e3, 90e3, 10e3, 64, line));
    V3 axis = vunit((V3){-m.g[1][1], -m.h[1][1], -m.g[1][0]});
    double err = 0;
    for (int k = 0; k <= 64; k++) {
        double z = vdot(vunit(line[k]), axis);
        err = fmax(err, fabs(vnorm(line[k]) / (1 - z * z) / (6371200. + 900e3) - 1));
    }
    CHECK(err < 3e-6);
    printf("dipole analytic field-line error %.8g\n", err);
    TubeGrid g;
    GridConfig c = {.nx = 4,
                    .ny = 4,
                    .ns = 32,
                    .lon = 114.2,
                    .span = 4,
                    .hmin = 250e3,
                    .hmax = 1400e3,
                    .hbase = 90e3,
                    .trace_step = 10e3};
    CHECK(!tube_geometry(&g, &m, c));
    CHECK(g.min_volume > 0 && g.max_parallel_electric_metric_error < 1e-13);
    tube_conductivity(&g, 1e10);
    double low = 1e100, high = 0;
    for (int i = 0; i < c.nx * c.ny * c.ns; i++) {
        CHECK(g.cell[i].pE >= 0 && g.cell[i].pF > 0);
        low = fmin(low, g.cell[i].pE);
        high = fmax(high, g.cell[i].pE);
    }
    CHECK(high > 0);
    printf("geometry positive volumes, E|| metric %.3g, field alignment %.3g\n",
           g.max_parallel_electric_metric_error, g.max_alignment_error);
    double *vacuum = calloc(c.nx * c.ny * c.ns, sizeof(double)),
           *coeff = calloc(6 * c.nx * c.ny, sizeof(double));
    potential_coefficients(&g, vacuum, 0, c.nx * c.ny, 1, coeff);
    for (int i = 0; i < c.nx * c.ny; i++) {
        double *k = coeff + 6 * i, symmetric = .5 * (k[1] + k[2]);
        CHECK(k[0] > 0 && k[3] > 0 && k[0] * k[3] > symmetric * symmetric);
    }
    puts(
        "F-fluid vacuum: integrated E-region Pedersen tensor remains positive without a conductivity floor");
    free(vacuum);
    free(coeff);
    tube_free(&g);
    double e1 = mms(12), e2 = mms(24), e3 = mms(48);
    CHECK(e1 / e2 > 3.6 && e2 / e3 > 3.6);
    regional_mms(12);
    e1 = regional_mms(24);
    e2 = regional_mms(48);
    e3 = regional_mms(96);
    CHECK(e1 / e2 > 3.4 && e2 / e3 > 3.4);
    e1 = periodic_mms(12);
    e2 = periodic_mms(24);
    e3 = periodic_mms(48);
    CHECK(e1 / e2 > 3.4 && e2 / e3 > 3.4);
    CHECK(!igrf_load(&m, argv[1], 2025, 13));
    double s64 = E_integral(&m, 64), s128 = E_integral(&m, 128), s512 = E_integral(&m, 512);
    printf("physical E-line Pedersen integral ns=64,128,512: %.12g %.12g %.12g S\n", s64, s128,
           s512);
    CHECK(fabs(s128 - s512) < .4 * fabs(s64 - s512));
    puts("CORE TESTS PASSED");
    return 0;
}
