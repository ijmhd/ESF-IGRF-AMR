#ifndef ESF_NEUTRAL_GRID_H
#define ESF_NEUTRAL_GRID_H
#include <stdint.h>
#include <stddef.h>
#include "tube.h"
/* Little-endian ESFNEU01, bound to the ordered canonical cell positions.
 * Per cell: O,N2,O2 [m^-3], Tn [K], nu_O,nu_e [s^-1]. Frozen bath.
 * Ti,Te in the header specify the constant temperatures used for collisions.
 * Return 0 for a legacy CSV, 1 on success, -1 for invalid binary/input. */
int neutral_grid_read(const char *path, const TubeGrid *g, double Ti, double Te, double **values);
uint64_t neutral_grid_hash(uint64_t h, const void *p, size_t n);
void tube_conductivity_values(TubeGrid *g, double e_peak, const double *values);
#endif
