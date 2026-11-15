#include "electric_driver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "electric driver check line%d\n", __LINE__);                           \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
/**
 * Check CSV validation, interpolation, coverage, knot clocks and driver identity.
 */
int main(int argc, char **argv) {
    CHECK(argc == 2);
    FILE *f = fopen(argv[1], "w");
    CHECK(f);
    fputs("time_s,east_V_m\n0,0\n3.7,0.0004\n9.2,-0.0002\n17.5,0\n", f);
    fclose(f);
    ElectricDriver *d = calloc(1, sizeof(*d));
    CHECK(d);
    CHECK(!electric_driver_load(d, argv[1], 17.5));
    CHECK(d->n == 4 && d->hash);
    CHECK(electric_driver_value(d, 3.7) == .0004);
    CHECK(fabs(electric_driver_value(d, 1.85) - .0002) < 1e-19);
    CHECK(fabs(electric_driver_value(d, 6.45) - .0001) < 1e-18);
    CHECK(electric_driver_next(d, 3.7) == 9.2);
    CHECK(isinf(electric_driver_next(d, 17.5)));
    CHECK(isnan(electric_driver_value(d, -1)));
    CHECK(isnan(electric_driver_value(d, 18)));
    CHECK(electric_driver_load(d, argv[1], 20));
    f = fopen(argv[1], "w");
    CHECK(f);
    fputs("time_s,east_V_m\n0,0\n0,1\n", f);
    fclose(f);
    CHECK(electric_driver_load(d, argv[1], 1));
    remove(argv[1]);
    free(d);
    puts("PASS driver interpolation, reversals, knot clocks, coverage and invalid input rejection");
    return 0;
}
