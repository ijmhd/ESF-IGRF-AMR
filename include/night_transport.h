#ifndef ESF_NIGHT_TRANSPORT_H
#define ESF_NIGHT_TRANSPORT_H
/* Isothermal, singly charged ions. SI arguments; no density/conductivity floor. */
void night_mobility(double nu, double omega, double *hall, double *pedersen);
double night_log_jump(double left, double right);
void night_sg_coeff(double Ti, double Te, double mass, double nu, double gravity, double ds,
                    double ne_left, double ne_right, double *theta, double *left, double *right);
double night_pressure_force(double Ti, double Te, double mass, double grad_i, double grad_e);
#endif
