/*
 * Local species mobility, pressure and collisional parallel-flux coefficients. Exact vacuum has
 * a bounded secant continuation rather than a positive density floor.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "night_transport.h"
#include <math.h>
static const double kb = 1.380649e-23;
/**
 * Evaluate bounded Hall/Pedersen mobility factors with hypot-based scaling.
 */
void night_mobility(double nu, double omega, double *hall, double *pedersen) {
    double h = hypot(nu, omega);
    *hall = (omega / h) * (omega / h);
    *pedersen = (nu / h) * (omega / h);
}
/**
 * Return log(right)-log(left) for positive states and a bounded secant continuation only at
 * exact vacuum.
 */
double night_log_jump(double left, double right) {
    if (left > 0 && right > 0)
        return log(right) - log(left);
    /* The pressure force at an exactly empty cell has no logarithm. Use the
     * bounded secant extension only at exact vacuum, never a density floor. */
    return left + right > 0 ? (right - left) / fmax(left, right) : 0;
}
/**
 * Evaluate the fitted-flux Bernoulli function robustly near zero and at large positive/negative
 * arguments.
 */
static double B(double u) {
    if (fabs(u) < 1e-4) {
        double u2 = u * u;
        return 1 - u / 2 + u2 / 12 - u2 * u2 / 720;
    }
    if (u > 50)
        return u * exp(-u) / (1 - exp(-u));
    if (u < -50)
        return -u / (1 - exp(u));
    return u / expm1(u);
}
/**
 * Compute collisional relaxation and Scharfetter-Gummel-style parallel pressure/gravity
 * coefficients from current electron density.
 */
void night_sg_coeff(double Ti, double Te, double mass, double nu, double gravity, double ds,
                    double ne_left, double ne_right, double *theta, double *left, double *right) {
    *theta = 1;
    *left = *right = 0;
    if (nu == 0)
        return;
    double ci2 = kb * Ti / mass, cs = sqrt(kb * (Ti + Te) / mass), kn = nu * ds / cs;
    *theta = 1 / (1 + kn * kn);
    double u = Te / Ti * night_log_jump(ne_left, ne_right) - gravity * ds / ci2;
    double d = (1 - *theta) * ci2 / (nu * ds);
    *left = d * B(u);
    *right = d * B(-u);
}
/**
 * Return the parallel pressure acceleration kB/m*(Ti*grad(log ni)+Te*grad(log Ne)).
 */
double night_pressure_force(double Ti, double Te, double mass, double grad_i, double grad_e) {
    return kb / mass * (Ti * grad_i + Te * grad_e);
}
