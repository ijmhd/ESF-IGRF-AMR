/*
 * Small local chemistry diagnostic used to inspect rate/charge behavior independently of spatial
 * transport.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "night.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/**
 * Parse this executable's interface, validate inputs, run its diagnostic/model entry point, and
 * return a nonzero status on failure.
 */
int main(int argc, char **argv) {
    if (argc != 10)
        return 2;
    double x[9];
    for (int i = 0; i < 9; i++) {
        char *end;
        x[i] = strtod(argv[i + 1], &end);
        if (*end || !isfinite(x[i]) || x[i] < 0)
            return 2;
    }
    if (x[8] <= 0)
        return 2;
    double o = x[0], no = x[1], o2 = x[2], t = 0, converted = 0, recombined = 0;
    NightRates r = {x[3], x[4], x[5], x[6]};
    while (t < x[7]) {
        double dt = fmin(x[8], x[7] - t);
        NightBudget b;
        if (night_react(&o, &no, &o2, r, dt, &b))
            return 1;
        converted += b.converted;
        recombined += b.recombined;
        t += dt;
    }
    printf("{\"O\":%.17g,\"NO\":%.17g,\"O2\":%.17g,\"converted\":%.17g,\"recombined\":%.17g}\n", o,
           no, o2, converted, recombined);
    return 0;
}
