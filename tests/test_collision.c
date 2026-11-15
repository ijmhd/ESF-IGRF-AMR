#include "tube.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "collision profile check failed line %d\n", __LINE__);                 \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)
/**
 * Check the prescribed collision/background response and input validation.
 */
int main(void) {
    const char *path = "collision_test.csv";
    FILE *f = fopen(path, "w");
    CHECK(f);
    fputs("height_m,nu_O_s1,nu_e_s1\n250000,2,200\n350000,.02,2\n", f);
    fclose(f);
    TubeCell c = {.height = 300e3,
                  .position = {6671200, 0, 0},
                  .B = 3e-5,
                  .b = {0, 0, 1},
                  .gx = {1, 0, 0},
                  .gy = {0, 1, 0}};
    TubeGrid g = {.c = {.nx = 1, .ny = 1, .ns = 1}, .cell = &c};
    CHECK(!tube_conductivity_table(&g, 0, path));
    double qe = 1.602176634e-19, mi = 16 * 1.66053906660e-27, me = 9.1093837139e-31;
    double ki = qe * c.B / mi / .2, ke = qe * c.B / me / 20;
    double expected = qe / c.B * 1e12 * (ki / (1 + ki * ki) + ke / (1 + ke * ke));
    CHECK(fabs(c.nu / .2 - 1) < 1e-14 && fabs(c.pF / expected - 1) < 1e-14 && c.pE == 0);
    double grav = -9.8 * pow(6371200 / 6671200., 2), oi = qe * c.B / mi, oe = qe * c.B / me;
    CHECK(fabs(c.cF[0] / (qe * 1e12 * grav * (.2 / (.2 * .2 + oi * oi) - 20 / (400 + oe * oe))) -
               1) < 1e-13);
    c.height = 240e3;
    CHECK(tube_conductivity_table(&g, 0, path) == -3);
    f = fopen(path, "w");
    CHECK(f);
    fputs("250000,2,200\n250000,.02,2\n", f);
    fclose(f);
    CHECK(tube_conductivity_table(&g, 0, path) == -2);
    remove(path);
    puts(
        "PASS collision input: consistent drag, conductivity and gravity current; coverage and ordering checks");
    return 0;
}
