#include "night_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "night transport check failed line %d\n", __LINE__);                   \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
/**
 * Check mobility, common-electron pressure and fitted parallel-flux limits/equilibria.
 */
int main(void) {
    double h, p;
    night_mobility(0, 100, &h, &p);
    CHECK(h == 1 && p == 0);
    night_mobility(1e12, 100, &h, &p);
    CHECK(h < 1.01e-20 && fabs(p / 1e-10 - 1) < 1e-12);
    night_mobility(100, 100, &h, &p);
    CHECK(fabs(h - .5) < 1e-15 && fabs(p - .5) < 1e-15);
    double Ti = 1000, Te = 1500, m = 30 * 1.66053906660e-27, g = -8, ds = 1000, th, l, r;
    /* Exact discrete hydrostatic equilibrium for a single ion Ne=ni. */
    double jump = m * g * ds / (1.380649e-23 * (Ti + Te));
    night_sg_coeff(Ti, Te, m, 1e3, g, ds, 1, exp(jump), &th, &l, &r);
    CHECK(fabs(l - r * exp(jump)) < 1e-13 * fmax(l, r));
    /* A molecular reservoir changes electron pressure, independent of ni. */
    night_sg_coeff(Ti, Te, m, 1e3, 0, ds, 1, 4, &th, &l, &r);
    CHECK(fabs(l / r - pow(4, -Te / Ti)) < 1e-14);
    double force = night_pressure_force(Ti, Te, m, 0, 2e-5);
    CHECK(fabs(force / (1.380649e-23 * Te / m * 2e-5) - 1) < 1e-15);
    CHECK(fabs(night_pressure_force(Ti, Te, m, jump / ds, jump / ds) - g) < 1e-12);
    night_sg_coeff(Ti, Te, m, 0, g, ds, 1, 2, &th, &l, &r);
    CHECK(th == 1 && l == 0 && r == 0);
    CHECK(isfinite(night_log_jump(0, 0)) && night_log_jump(0, 0) == 0);
    puts(
        "PASS species mobility limits, common electron pressure, SG hydrostatic equilibrium and vacuum handling");
}
