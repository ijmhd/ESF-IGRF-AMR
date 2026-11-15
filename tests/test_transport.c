#include "solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PI 3.14159265358979323846
#define CHECK(c)                                                                                   \
    do {                                                                                           \
        if (!(c)) {                                                                                \
            fprintf(stderr, "transport check failed line %d\n", __LINE__);                         \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
/**
 * Build a synthetic Cartesian finite-volume geometry for transport tests.
 */
static void box(Solver3 *s, int n, int nz, double length) {
    s->c.grid.nx = s->c.grid.ny = n;
    s->c.grid.ns = nz;
    s->g.c = s->c.grid;
    s->size = 1;
    s->count = n * n;
    s->n = s->total = n * n * nz;
    s->g.cell = calloc(s->total, sizeof(TubeCell));
    s->rho = calloc(s->total, sizeof(double));
    s->momentum = calloc(s->total, sizeof(double));
    s->stage_n = calloc(s->total, sizeof(double));
    s->stage_m = calloc(s->total, sizeof(double));
    s->old_n = calloc(s->total, sizeof(double));
    s->old_m = calloc(s->total, sizeof(double));
    s->rhs_n = calloc(s->total, sizeof(double));
    s->rhs_m = calloc(s->total, sizeof(double));
    s->phi = calloc((n + 1) * (n + 1), sizeof(double));
    s->first_rhs_m = calloc(s->total, sizeof(double));
    for (int i = 0; i < s->total; i++)
        s->g.cell[i].volume = length / s->total;
    for (int a = 0; a < 3; a++) {
        int count = tube_face_count(&s->g, a);
        s->g.face[a] = calloc(count, sizeof(TubeFace));
        s->flux_n[a] = calloc(count, sizeof(double));
        s->flux_m[a] = calloc(count, sizeof(double));
        s->rate_face[a] = calloc(count, sizeof(double));
        for (int k = 0; k < count; k++) {
            if (a == 0)
                s->g.face[a][k].e[1] = -length / (n * nz);
            if (a == 1)
                s->g.face[a][k].e[0] = length / (n * nz);
            if (a == 2)
                s->g.face[a][k].parallel_area = 1. / (n * n);
        }
    }
}
/**
 * Check constant-state behavior of the uniform transport operator.
 */
static void flat(Solver3 *s, int n) {
    box(s, n, 4, 1);
}
/**
 * Check the selected small-amplitude parallel acoustic fixture.
 */
static double acoustic(int nz, int order, double nu, int relaxation, int prepared) {
    Solver3 s = {0};
    const double L = 1e6, cs = 1000, eps = 1e-5, k = PI / L, end = 1000;
    box(&s, 4, nz, L);
    s.c.parallel = 1;
    s.c.parallel_order = order;
    s.c.parallel_relaxation = relaxation;
    s.c.temperature = 16 * 1.66053906660e-27 * cs * cs / 1.380649e-23;
    for (int p = 0; p < 16; p++)
        for (int z = 0; z < nz; z++) {
            int id = p * nz + z;
            double phase = k * (z + .5) * L / nz,
                   average = sin(.5 * k * L / nz) / (.5 * k * L / nz);
            s.rho[id] = 1 + eps * average * cos(phase);
            s.g.cell[id].nu = nu;
            if (nu && prepared) {
                double root = sqrt(nu * nu - 4 * cs * cs * k * k),
                       slow = -2 * cs * cs * k * k / (nu + root);
                s.momentum[id] = -eps * slow / k * average * sin(phase);
            }
        }
    double time = 0, mass = solver3_mass(&s, s.rho);
    while (time < end - 1e-10) {
        solver3_transport(&s, s.rho, s.momentum, s.rhs_n, s.rhs_m);
        double dt = fmin(.3 / s.cfl_rate, end - time);
        for (int id = 0; id < s.total; id++) {
            double d, p1, p2;
            exponential_weights(nu, dt, &d, &p1, &p2);
            s.old_n[id] = s.rho[id];
            s.old_m[id] = s.momentum[id];
            s.first_rhs_m[id] = s.rhs_m[id];
            s.stage_n[id] = s.rho[id] + dt * s.rhs_n[id];
            s.stage_m[id] = d * s.momentum[id] + dt * p1 * s.rhs_m[id];
        }
        solver3_transport(&s, s.stage_n, s.stage_m, s.rhs_n, s.rhs_m);
        for (int id = 0; id < s.total; id++) {
            double d, p1, p2;
            exponential_weights(nu, dt, &d, &p1, &p2);
            s.rho[id] = .5 * s.old_n[id] + .5 * (s.stage_n[id] + dt * s.rhs_n[id]);
            s.momentum[id] =
                d * s.old_m[id] + dt * ((p1 - p2) * s.first_rhs_m[id] + p2 * s.rhs_m[id]);
            CHECK(isfinite(s.rho[id]) && s.rho[id] > 0);
        }
        time += dt;
    }
    double factor = cos(cs * k * end);
    if (nu) {
        double root = sqrt(nu * nu - 4 * cs * cs * k * k),
               slow = -2 * cs * cs * k * k / (nu + root), fast = -.5 * (nu + root);
        factor = prepared ? exp(slow * end)
                          : (-fast * exp(slow * end) + slow * exp(fast * end)) / (slow - fast);
    }
    double err = 0;
    for (int z = 0; z < nz; z++) {
        double average = sin(.5 * k * L / nz) / (.5 * k * L / nz),
               truth = 1 + eps * average * cos(k * (z + .5) * L / nz) * factor;
        err += fabs(s.rho[z] - truth) / (eps * nz);
    }
    CHECK(fabs(solver3_mass(&s, s.rho) - mass) / mass < 2e-13);
    printf(
        "linear isothermal acoustic mode ns=%d order=%d nu=%.4g relaxation=%d prepared=%d perturbation_L1=%.9g exact_amplitude=%.9g\n",
        nz, order, nu, relaxation, prepared, err, factor);
    solver3_free(&s);
    return err;
}
/**
 * Check the idealized initialization/transport state against its declared fixture.
 */
static double initial(double x, double y, double z) {
    double d = hypot(x - .5, y - .60) / .10;
    return 1 + (d < 1 ? .1 * pow(1 - d * d, 4) * (1 + .2 * cos(2 * PI * z)) : 0);
}
/**
 * Check the rotational transport fixture and its conservation diagnostics.
 */
static double rotation(int n, int pdm) {
    Solver3 s = {0};
    flat(&s, n);
    const double q[4] = {-.8611363115940526, -.3399810435848563, .3399810435848563,
                         .8611363115940526},
                 w[4] = {.3478548451374538, .6521451548625461, .6521451548625461,
                         .3478548451374538};
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            for (int k = 0; k < 4; k++) {
                double value = 0;
                for (int a = 0; a < 4; a++)
                    for (int b = 0; b < 4; b++)
                        for (int c = 0; c < 4; c++)
                            value += .125 * w[a] * w[b] * w[c] *
                                     initial((i + .5 + .5 * q[a]) / n, (j + .5 + .5 * q[b]) / n,
                                             (k + .5 + .5 * q[c]) / 4);
                int id = (j * n + i) * 4 + k;
                s.rho[id] = s.old_n[id] = value;
            }
    if (pdm) {
        s.c.exb_pdm4 = 1;
        s.c.pdm_beta = 1;
        s.g.node = calloc((size_t)(n + 1) * (n + 1) * 5, sizeof(V3));
        for (int j = 0; j <= n; j++)
            for (int i = 0; i <= n; i++)
                for (int k = 0; k <= 4; k++)
                    s.g.node[(j * (n + 1) + i) * 5 + k] =
                        (V3){(double)i / n, (double)j / n, (double)k / 4};
        for (int a = 0; a < 3; a++)
            s.pdm_weights[a] = calloc(4 * tube_face_count(&s.g, a), sizeof(double));
        CHECK(!tube_pdm_weights(&s.g, 0, s.pdm_weights));
    }
    /* Quadratic streamfunction in the compact blob's orbit; C1 taper to
     * constant outside r=.4. The perturbation stays entirely inside r=.25. */
    for (int j = 0; j <= n; j++)
        for (int i = 0; i <= n; i++) {
            double r = hypot((double)i / n - .5, (double)j / n - .5), f;
            if (r <= .25)
                f = .5 * r * r;
            else if (r >= .4)
                f = .5 * .25 * .25 + .25 * .15 / 2;
            else {
                double t = (r - .25) / .15;
                f = .5 * .25 * .25 + .25 * .15 * (t - t * t / 2);
            }
            s.phi[j * (n + 1) + i] = 2 * PI * f;
        }
    /* A constant scalar must remain constant under the discrete curl. */
    double *constant = calloc(s.total, sizeof(double));
    for (int k = 0; k < s.total; k++)
        constant[k] = 1;
    CHECK(!solver3_transport(&s, constant, s.momentum, s.rhs_n, s.rhs_m));
    double err = 0;
    for (int k = 0; k < s.total; k++)
        err = fmax(err, fabs(s.rhs_n[k]));
    CHECK(err < 1e-11);
    free(constant);
    double time = 0, mass = solver3_mass(&s, s.rho);
    const double end = .125;
    while (time < end - 1e-13) {
        solver3_transport(&s, s.rho, s.momentum, s.rhs_n, s.rhs_m);
        double dt = fmin(.35 / s.cfl_rate, end - time);
        for (int k = 0; k < s.total; k++)
            s.stage_n[k] = s.rho[k] + dt * s.rhs_n[k];
        solver3_transport(&s, s.stage_n, s.momentum, s.rhs_n, s.rhs_m);
        for (int k = 0; k < s.total; k++)
            s.rho[k] = .5 * s.rho[k] + .5 * (s.stage_n[k] + dt * s.rhs_n[k]);
        time += dt;
    }
    err = 0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            for (int k = 0; k < 4; k++) {
                double truth = 0;
                for (int a = 0; a < 4; a++)
                    for (int b = 0; b < 4; b++)
                        for (int c = 0; c < 4; c++) {
                            double x = (i + .5 + .5 * q[a]) / n - .5,
                                   y = (j + .5 + .5 * q[b]) / n - .5, angle = 2 * PI * end;
                            truth += .125 * w[a] * w[b] * w[c] *
                                     initial(.5 + cos(angle) * x + sin(angle) * y,
                                             .5 - sin(angle) * x + cos(angle) * y,
                                             (k + .5 + .5 * q[c]) / 4);
                        }
                err += fabs(s.rho[(j * n + i) * 4 + k] - truth) / s.total;
            }
    CHECK(fabs(solver3_mass(&s, s.rho) - mass) < 2e-13);
    printf("3D-volume scalar rotation %d pdm=%d error %.12g\n", n, pdm, err);
    solver3_free(&s);
    return err;
}
/**
 * Check uniform finite-volume constant, acoustic, initial-profile and rotational transport
 * fixtures.
 */
int main(void) {
    Solver3 s = {0};
    flat(&s, 8);
    s.c.parallel = 1;
    s.c.temperature = 2500;
    for (int i = 0; i < s.total; i++)
        s.rho[i] = 1.;
    for (int j = 0; j < 8; j++)
        for (int i = 0; i < 8; i++)
            for (int k = 0; k <= 4; k++)
                s.g.face[2][tube_face_id(&s.g, 2, i, j, k)].parallel_area *=
                    (1 + .2 * sin(2 * PI * k / 4));
    CHECK(!solver3_transport(&s, s.rho, s.momentum, s.rhs_n, s.rhs_m));
    double e = 0;
    for (int i = 0; i < s.total; i++) {
        CHECK(s.rhs_n[i] == 0);
        e = fmax(e, fabs(s.rhs_m[i]));
    }
    CHECK(e < 1e-8);
    printf("variable-area isothermal rest-state acceleration %.6g\n", e);
    solver3_free(&s);
    double a = rotation(32, 0), b = rotation(64, 0);
    CHECK(a / b > 2.5);
    CHECK(rotation(32, 1) < a);
    a = acoustic(32, 2, 0, 0, 0);
    b = acoustic(64, 2, 0, 0, 0);
    CHECK(a / b > 2.8);
    a = acoustic(32, 2, 0, 1, 0);
    b = acoustic(64, 2, 0, 1, 0);
    CHECK(a / b > 2.8);
    a = acoustic(64, 1, 1, 0, 0);
    b = acoustic(64, 2, 1, 0, 0);
    CHECK(b < .4 * a);
    a = acoustic(64, 2, 1000, 0, 0);
    b = acoustic(64, 2, 1000, 1, 0);
    CHECK(b < 1e-7 && b < .001 * a);
    a = acoustic(32, 2, .01, 1, 0);
    b = acoustic(64, 2, .01, 1, 0);
    CHECK(a / b > 2.5);
    /* The under-resolved initial layer has an O(1/nu^2) error plateau; keep
     * its independent absolute check, and test space convergence separately
     * on a well-prepared slow eigenmode instead of hiding the plateau. */
    CHECK(acoustic(64, 2, 1, 1, 0) < 1e-5);
    a = acoustic(32, 2, 1, 1, 1);
    b = acoustic(64, 2, 1, 1, 1);
    CHECK(a / b > 2.5);
    Solver3 eq = {0};
    box(&eq, 4, 64, 1e6);
    eq.c.parallel = 1;
    eq.c.parallel_order = 2;
    eq.c.parallel_relaxation = 1;
    eq.c.temperature = 16 * 1.66053906660e-27 * 1e6 / 1.380649e-23;
    for (int p = 0; p < 16; p++)
        for (int k = 0; k < 64; k++) {
            int id = p * 64 + k;
            double x = (k + .5) / 64.;
            eq.rho[id] = exp(x) * sinh(.5 / 64) / (.5 / 64);
            eq.g.cell[id].nu = 1e6;
            eq.g.cell[id].gravity_parallel = 1;
        }
    solver3_transport(&eq, eq.rho, eq.momentum, eq.rhs_n, eq.rhs_m);
    double en = 0, em = 0;
    for (int id = 0; id < eq.total; id++) {
        en = fmax(en, fabs(eq.rhs_n[id]));
        em = fmax(em, fabs(eq.rhs_m[id]));
    }
    printf(
        "isothermal exponential hydrostatic state, strong drag: density RHS %.9g momentum RHS %.9g\n",
        en, em);
    CHECK(en < 1e-12 && em < 1e-10);
    solver3_free(&eq);
    puts("TRANSPORT TESTS PASSED");
    return 0;
}
