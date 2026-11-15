/*
 * Strict reader for externally supplied frozen neutral/collision grids. Runtime does not
 * evaluate MSIS or an energy equation.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "neutral_grid.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/**
 * Compute the frozen-neutral input file checksum used in physical restart/cache signatures.
 */
uint64_t neutral_grid_hash(uint64_t h, const void *p, size_t n) {
    const unsigned char *b = p;
    while (n--) {
        h ^= *b++;
        h *= UINT64_C(1099511628211);
    }
    return h;
}
/**
 * Validate binary magic, dimensions, positions, temperatures, finite neutral/collision data and
 * checksum against canonical geometry; return owned allocated input values.
 */
int neutral_grid_read(const char *path, const TubeGrid *g, double Ti, double Te, double **values) {
    *values = NULL;
    FILE *f = fopen(path, "rb");
    if (!f)
        return -1;
    char magic[8];
    if (fread(magic, 1, 8, f) != 8) {
        fclose(f);
        return -1;
    }
    if (memcmp(magic, "ESFNEU01", 8)) {
        fclose(f);
        return 0;
    }
    unsigned short endian = 1;
    uint64_t h[7];
    double temperatures[2];
    if (!*(unsigned char *)&endian || sizeof(double) != 8 || fread(h, 8, 7, f) != 7 ||
        fread(temperatures, 8, 2, f) != 2 || h[0] != 1 || h[1] != (uint64_t)g->c.nx ||
        h[2] != (uint64_t)g->c.ny || h[3] != (uint64_t)g->c.ns || h[4] != 6 ||
        temperatures[0] != Ti || temperatures[1] != Te) {
        fclose(f);
        return -1;
    }
    size_t n = (size_t)g->c.nx * g->c.ny * g->c.ns;
    if (n > SIZE_MAX / (6 * sizeof(double))) {
        fclose(f);
        return -1;
    }
    uint64_t position = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < n; i++)
        position = neutral_grid_hash(position, &g->cell[i].position, 24);
    if (position != h[5]) {
        fclose(f);
        return -1;
    }
    double *v = malloc(6 * n * sizeof(double));
    if (!v) {
        fclose(f);
        return -1;
    }
    int good = fread(v, sizeof(double), 6 * n, f) == 6 * n && fgetc(f) == EOF && !ferror(f);
    if (fclose(f))
        good = 0;
    if (good)
        good = neutral_grid_hash(UINT64_C(14695981039346656037), v, 6 * n * sizeof(double)) == h[6];
    for (size_t i = 0; good && i < 6 * n; i++)
        if (!isfinite(v[i]) || v[i] <= 0)
            good = 0;
    if (!good) {
        free(v);
        return -1;
    }
    *values = v;
    return 1;
}
