#ifndef ESF_INTEGRATOR_H
#define ESF_INTEGRATOR_H
/* ETD2 weights for y' = -nu*y + R(y), nu fixed and nonnegative.
 * The exact drag/constant-force balance is retained even for nu*dt >> 1. */
void exponential_weights(double nu, double dt, double *decay, double *phi1, double *phi2);
#endif
