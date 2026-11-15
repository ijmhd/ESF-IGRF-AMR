/*
 * Validated piecewise-linear imposed eastward electric field. Times are elapsed model seconds;
 * field units in the CSV are V/m.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "electric_driver.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
/**
 * Read/hash a strict time_s,east_V_m CSV, require finite increasing times from zero and coverage
 * through the requested absolute end time.
 */
int electric_driver_load(ElectricDriver *d, const char *path, double end) {
    memset(d, 0, sizeof(*d));
    d->hash = UINT64_C(14695981039346656037);
    FILE *f = fopen(path, "rb");
    if (!f)
        return 1;
    int ch;
    while ((ch = fgetc(f)) != EOF) {
        d->hash ^= (unsigned char)ch;
        d->hash *= UINT64_C(1099511628211);
    }
    if (ferror(f)) {
        fclose(f);
        return 2;
    }
    rewind(f);
    char line[256], extra;
    if (!fgets(line, sizeof(line), f) || strcmp(line, "time_s,east_V_m\n")) {
        fclose(f);
        return 3;
    }
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n')
            continue;
        int n = d->n;
        if (n == ELECTRIC_DRIVER_MAX ||
            sscanf(line, "%lf,%lf %c", d->t + n, d->east + n, &extra) != 2 || !isfinite(d->t[n]) ||
            !isfinite(d->east[n]) || d->t[n] < 0 || (n && d->t[n] <= d->t[n - 1])) {
            fclose(f);
            return 4;
        }
        d->n++;
    }
    int error = ferror(f) || d->n < 2 || d->t[0] != 0 || d->t[d->n - 1] < end;
    fclose(f);
    return error ? 5 : 0;
}
/**
 * Linearly interpolate the imposed field at a model stage time; return NaN outside declared
 * coverage.
 */
double electric_driver_value(const ElectricDriver *d, double time) {
    if (time < d->t[0] - 1e-9 || time > d->t[d->n - 1] + 1e-9)
        return NAN;
    if (time <= d->t[0])
        return d->east[0];
    if (time >= d->t[d->n - 1])
        return d->east[d->n - 1];
    int lo = 0, hi = d->n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (d->t[mid] > time)
            hi = mid;
        else
            lo = mid;
    }
    double w = (time - d->t[lo]) / (d->t[hi] - d->t[lo]);
    return (1 - w) * d->east[lo] + w * d->east[hi];
}
/**
 * Return the next forcing knot after time, or infinity if absent, so accepted steps land on
 * forcing changes.
 */
double electric_driver_next(const ElectricDriver *d, double time) {
    if (d)
        for (int i = 0; i < d->n; i++)
            if (d->t[i] > time + 1e-9)
                return d->t[i];
    return INFINITY;
}
