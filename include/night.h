#ifndef ESF_NIGHT_H
#define ESF_NIGHT_H
/* SI densities/rates. Four-reaction nighttime O+/NO+/O2+ kernel.
 * Rate fits: GEMINI3D sources.f90, commit d70446a712d47e7617b25c033bfcace64207590f.
 * Thermal rates assume zero ion-neutral relative drift; not full photochemistry. */
typedef struct {
    double to_no, to_o2, alpha_no, alpha_o2;
} NightRates;
typedef struct {
    double converted, recombined;
} NightBudget;
int night_rates(double Ti, double Te, double Tn, double N2, double O2, NightRates *r);
int night_react(double *o, double *no, double *o2, NightRates r, double dt, NightBudget *budget);
double night_solar_cosine(double unix_utc, double latitude_rad, double longitude_rad);
#endif
