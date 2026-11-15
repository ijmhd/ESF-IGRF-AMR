#ifndef ESF_ELECTRIC_DRIVER_H
#define ESF_ELECTRIC_DRIVER_H
#include <stdint.h>
#define ELECTRIC_DRIVER_MAX 4096
typedef struct {
    int n;
    double t[ELECTRIC_DRIVER_MAX], east[ELECTRIC_DRIVER_MAX];
    uint64_t hash;
} ElectricDriver;
/* CSV time_s,east_V_m; SI electric field at the declared apex reference. */
int electric_driver_load(ElectricDriver *d, const char *path, double end);
double electric_driver_value(const ElectricDriver *d, double time);
double electric_driver_next(const ElectricDriver *d, double time);
#endif
