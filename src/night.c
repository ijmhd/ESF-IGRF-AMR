/*
 * Local reaction and solar-geometry kernels; no MPI or mutable global state. Reaction density
 * units are m^-3, rates s^-1 or m^3 s^-1, and time seconds.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "night.h"
#include <math.h>
#define PI 3.14159265358979323846
/**
 * Calculate O+-to-NO+/O2+ conversion and dissociative recombination rates from prescribed
 * temperatures and N2/O2; reject unsupported thermal ranges.
 */
int night_rates(double Ti, double Te, double Tn, double N2, double O2, NightRates *r) {
    if (!r || !isfinite(Ti + Te + Tn + N2 + O2) || Ti < 200 || Ti > 5000 || Te < 200 || Te > 5000 ||
        Tn < 100 || Tn > 5000 || N2 < 0 || O2 < 0)
        return -1;
    double x = (28. / 44 * (Ti - Tn) + Tn) / 300.;
    /* Low-Teff branches valid for this declared thermal temperature range. */
    if (x * 300 > 3725)
        return -1;
    double k1 =
        1.71676e-12 + x * (-7.19934e-13 + x * (1.33276e-13 + x * (-9.28213e-15 + x * 6.39557e-16)));
    x = (32. / 48 * (Ti - Tn) + Tn) / 300.;
    if (x * 300 > 4800)
        return -1;
    double k2 =
        2.78932e-11 + x * (-6.92612e-12 + x * (8.67684e-13 + x * (-3.47251e-14 + x * 5.07097e-16)));
    r->to_no = k1 * 1e-6 * N2;
    r->to_o2 = k2 * 1e-6 * O2;
    r->alpha_no = 4.2e-13 * pow(300 / Te, .85);
    r->alpha_o2 = Te <= 1200 ? 1.95e-13 * pow(300 / Te, .70) : 7.38e-14 * pow(1200 / Te, .56);
    return r->to_no >= 0 && r->to_o2 >= 0 ? 0 : -1;
}
/**
 * Apply the cached exact linear conversion operator, splitting lost O+ between molecular
 * products and accumulating converted particles.
 */
static void convert(double *o, double *no, double *o2, double fraction, double decay,
                    double branch_no, double branch_o2, NightBudget *b) {
    double loss = *o * fraction;
    *o *= decay;
    *no += loss * branch_no;
    *o2 += loss * branch_o2;
    b->converted += loss;
}
/**
 * Apply the exact scalar solution of dn/dt=-alpha*n*(other+n) for fixed other species, with a
 * stable zero-other limit.
 */
static double recombine(double n, double other, double alpha, double dt) {
    /* Exact dn/dt=-alpha*n*(other+n); stable also as other -> 0. */
    double x = alpha * other * dt;
    if (x == 0)
        return n / (1 + alpha * n * dt);
    double exprel = -expm1(-x) / x;
    return n * exp(-x) / (1 + alpha * n * dt * exprel);
}
/**
 * Compose symmetric conversion/NO+/O2+/NO+/conversion operators, retaining nonnegativity and an
 * explicit recombined-electron budget; no photoionization is included.
 */
int night_react(double *o, double *no, double *o2, NightRates r, double dt, NightBudget *b) {
    if (!o || !no || !o2 || !b || !isfinite(*o + *no + *o2 + dt) || *o < 0 || *no < 0 || *o2 < 0 ||
        dt < 0 || !isfinite(r.to_no + r.to_o2 + r.alpha_no + r.alpha_o2) || r.to_no < 0 ||
        r.to_o2 < 0 || r.alpha_no < 0 || r.alpha_o2 < 0)
        return -1;
    /* Both conversion half operators have the same rates and duration.
     * Reuse only these state-independent coefficients. Recombination below
     * depends on the updated electron density and must be reevaluated. */
    double rate = r.to_no + r.to_o2, fraction = 0, decay = 1, branch_no = 0, branch_o2 = 0;
    if (rate > 0) {
        double half = .5 * dt;
        fraction = -expm1(-rate * half);
        decay = exp(-rate * half);
        branch_no = r.to_no / rate;
        branch_o2 = r.to_o2 / rate;
    }
    *b = (NightBudget){0};
    if (rate > 0)
        convert(o, no, o2, fraction, decay, branch_no, branch_o2, b);
    double old = *no;
    *no = recombine(*no, *o + *o2, r.alpha_no, .5 * dt);
    b->recombined += old - *no;
    old = *o2;
    *o2 = recombine(*o2, *o + *no, r.alpha_o2, dt);
    b->recombined += old - *o2;
    old = *no;
    *no = recombine(*no, *o + *o2, r.alpha_no, .5 * dt);
    b->recombined += old - *no;
    if (rate > 0)
        convert(o, no, o2, fraction, decay, branch_no, branch_o2, b);
    return isfinite(*o + *no + *o2) && *o >= 0 && *no >= 0 && *o2 >= 0 ? 0 : -1;
}
/**
 * Estimate the solar zenith cosine from elapsed UTC and geographic angles for the declared
 * night-only guard, not precise sunrise or MLT.
 */
double night_solar_cosine(double utc, double lat, double lon) {
    /* Low-order solar ephemeris, UTC->Julian date, apparent solar direction.
     * Used for a night-only guard, not precision sunrise timing or MLT. */
    double jd = utc / 86400 + 2440587.5, d = jd - 2451545.;
    double L = fmod(280.460 + .9856474 * d, 360) * PI / 180,
           g = fmod(357.528 + .9856003 * d, 360) * PI / 180;
    double lambda = L + (1.915 * sin(g) + .020 * sin(2 * g)) * PI / 180,
           eps = (23.439 - .0000004 * d) * PI / 180;
    double ra = atan2(cos(eps) * sin(lambda), cos(lambda)), dec = asin(sin(eps) * sin(lambda));
    double gmst = (280.46061837 + 360.98564736629 * d) * PI / 180;
    return sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(gmst + lon - ra);
}
