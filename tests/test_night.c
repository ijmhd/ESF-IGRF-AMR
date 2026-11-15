#include "night.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define CHECK(x)                                                                                   \
    do {                                                                                           \
        if (!(x)) {                                                                                \
            fprintf(stderr, "night check failed line%d\n", __LINE__);                              \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
/**
 * Check local reactions, positivity, exact limiting behavior and charge accounting.
 */
int main(void) {
    NightRates r;
    CHECK(!night_rates(1000, 1500, 300, 1e16, 1e15, &r));
    double o = 1e11, no = 2e10, o2 = 3e9, total = o + no + o2;
    NightBudget b;
    CHECK(!night_react(&o, &no, &o2, (NightRates){.to_no = .01, .to_o2 = .02}, 100, &b));
    CHECK(fabs(o / (1e11 * exp(-3)) - 1) < 1e-14);
    CHECK(fabs((o + no + o2) / total - 1) < 1e-14);
    o = 0;
    no = 1e10;
    o2 = 0;
    double alpha = 4e-13;
    CHECK(!night_react(&o, &no, &o2, (NightRates){.alpha_no = alpha}, 2000, &b));
    CHECK(fabs(no / (1e10 / (1 + alpha * 1e10 * 2000)) - 1) < 1e-14);
    CHECK(fabs((1e10 - no - b.recombined) / 1e10) < 1e-14);
    o = 1e12;
    no = 1e10;
    o2 = 2e10;
    total = o + no + o2;
    CHECK(!night_react(&o, &no, &o2, (NightRates){1e4, 2e4, 1e-12, 2e-12}, 100, &b));
    CHECK(fabs((o + no + o2 + b.recombined) / total - 1) < 2e-14);
    CHECK(o >= 0 && no > 0 && o2 >= 0);
    CHECK(night_react(&o, &no, &o2, (NightRates){-1, 0, 0, 0}, 1, &b));
    /* Equinox near Greenwich noon: subsolar cosine~1, midnight~-1. */
    double utc = 1742472000.;
    CHECK(night_solar_cosine(utc, 0, 0) > .99);
    CHECK(night_solar_cosine(utc, 0, 3.141592653589793) < -.99);
    puts(
        "PASS night reaction positivity, conversion charge, recombination balance, stiff kinetics and solar clock");
    return 0;
}
