/*
 * Exponential time differencing for fixed nonnegative collision drag; density uses the matching
 * two-stage explicit update.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "integrator.h"
#include <math.h>
/**
 * Evaluate exp(-nu*dt), phi1 and phi2 using stable small-argument series and expm1. Reusing
 * these weights is valid only for the same nu and attempted dt.
 */
void exponential_weights(double nu, double dt, double *decay, double *p1, double *p2) {
    double a = nu * dt;
    if (a == 0) {
        *decay = 1;
        *p1 = 1;
        *p2 = .5;
        return;
    }
    *decay = exp(-a);
    if (a < 1e-4) {
        *p1 = 1 - a / 2 + a * a / 6 - a * a * a / 24 + a * a * a * a / 120;
        *p2 = .5 - a / 6 + a * a / 24 - a * a * a / 120 + a * a * a * a / 720;
    } else {
        *p1 = -expm1(-a) / a;
        *p2 = (1 - *p1) / a;
    }
}
