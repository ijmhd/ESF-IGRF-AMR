/*
 * Native-ABI geometry cache. This is a performance cache, not a portable state checkpoint.
 * Source/input identities, ABI sizes and payload checksum must all match.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Local native-ABI cache, never a portable checkpoint. Complete input/code
 * identity, ABI sizes and a payload hash are checked before accepting a cache.
 */
#include "amr3.h"
#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/**
 * Update the FNV-1a checksum over an exact byte sequence.
 */
static uint64_t hash(uint64_t h, const void *data, size_t n) {
    const unsigned char *p = data;
    while (n--) {
        h ^= *p++;
        h *= UINT64_C(1099511628211);
    }
    return h;
}
/**
 * Hash an input file for cache identity; an empty optional path has identity zero.
 */
static uint64_t file_hash(A3 *a, const char *path) {
    if (!path[0])
        return 0;
    FILE *f = fopen(path, "rb");
    if (!f)
        a3_fail(a, "cache input hash open");
    unsigned char block[8192];
    size_t n;
    uint64_t h = UINT64_C(14695981039346656037);
    while ((n = fread(block, 1, sizeof(block), f)))
        h = hash(h, block, n);
    if (ferror(f))
        a3_fail(a, "cache input read");
    fclose(f);
    return h;
}
/**
 * Accept the exact compiled geometry identity or one explicitly enumerated, version-guarded
 * historical equivalence; reject unknown changes.
 */
int a3_geometry_identity(A3 *a, const char code[128]) {
    if (!memcmp(code, ESF_GEOMETRY_ID, 128))
        return 1;
    /* The first parallel-geometry build retained the configure-time serial
     * source fingerprint. Only wrappers/row ranges changed, not arithmetic.
     * This exact pair was audited by old/new static and fine-grid equivalence.
     * Any future geometry/IGRF edit disables this compatibility automatically. */
    const char old[] =
        "ad6ed0de002e6be3f6bdaf83bdf947ca8acd0da64d3341ed246980dbbf035bcdc3f2071a51c2a9d6d9777937513eb9986c587718076689f420a1431d7df7e2db";
    const char equivalent[] =
        "880f1def6afb918e70abc212846c8042cd45046cfc9c07028443f11895937940c3f2071a51c2a9d6d9777937513eb9986c587718076689f420a1431d7df7e2db";
    if (!memcmp(ESF_GEOMETRY_ID, equivalent, 128) && !memcmp(code, old, 128)) {
        if (!a->rank) {
            fprintf(stdout, "Verified compatible legacy serial/parallel geometry identity\n");
            fflush(stdout);
        }
        return 1;
    }
    /* Only an appended conductivity API differs; tracing and metric arithmetic
     * are unchanged. Regression checks compare the complete old/new states. */
    const char before_spatial_api[] =
        "880f1def6afb918e70abc212846c8042cd45046cfc9c07028443f11895937940c3f2071a51c2a9d6d9777937513eb9986c587718076689f420a1431d7df7e2db";
    const char with_spatial_api[] =
        "c42f5f196152bb2fe7ac678ecb4eddb63cada4306ec146a09332a61b087d2c3dc3f2071a51c2a9d6d9777937513eb9986c587718076689f420a1431d7df7e2db";
    if (!memcmp(ESF_GEOMETRY_ID, with_spatial_api, 128) && !memcmp(code, before_spatial_api, 128)) {
        if (!a->rank) {
            fprintf(stdout, "Verified geometry identity for appended spatial conductivity API\n");
            fflush(stdout);
        }
        return 1;
    }
    return 0;
}
/**
 * Read/write native canonical geometry with complete parameter/input/source/ABI checks, atomic
 * rename and payload checksum.
 */
int a3_geometry_cache(A3 *a, TubeGrid *g, Config3 c, int write) {
    if (!c.geometry_cache[0])
        return 0;
    uint64_t id[] = {1,
                     sizeof(V3),
                     sizeof(TubeCell),
                     sizeof(TubeFace),
                     file_hash(a, c.coeff),
                     file_hash(a, c.collision_profile),
                     UINT64_C(0x0102030405060708)};
    GridConfig f = c.grid;
    f.nx *= 2;
    f.ny *= 2;
    double pars[] = {f.nx,
                     f.ny,
                     f.ns,
                     f.lon,
                     f.span,
                     f.hmin,
                     f.hmax,
                     f.hbase,
                     f.trace_step,
                     c.year,
                     c.degree,
                     c.aligned_dipole,
                     c.dipole_equator_nT,
                     c.e_peak};
    char temporary[1200];
    snprintf(temporary, sizeof(temporary), "%s.part.%ld", c.geometry_cache, (long)getpid());
    FILE *stream = fopen(write ? temporary : c.geometry_cache, write ? "wb" : "rb");
    if (!stream) {
        if (!write && errno == ENOENT)
            return 0;
        a3_fail(a, "geometry cache open");
    }
    if (write) {
        if (fwrite("ESFGEO01", 1, 8, stream) != 8 ||
            fwrite(id, 1, sizeof(id), stream) != sizeof(id) ||
            fwrite(pars, 1, sizeof(pars), stream) != sizeof(pars) ||
            fwrite(ESF_GEOMETRY_ID, 1, 128, stream) != 128)
            a3_fail(a, "cache header write");
    } else {
        char magic[8], code[128];
        uint64_t got[7];
        double params[14];
        if (fread(magic, 1, 8, stream) != 8 || memcmp(magic, "ESFGEO01", 8) ||
            fread(got, 1, sizeof(got), stream) != sizeof(got) || memcmp(got, id, sizeof(id)) ||
            fread(params, 1, sizeof(params), stream) != sizeof(params) ||
            memcmp(params, pars, sizeof(pars)) || fread(code, 1, 128, stream) != 128 ||
            !a3_geometry_identity(a, code))
            a3_fail(a, "geometry cache input/code/ABI mismatch");
        memset(g, 0, sizeof(*g));
        g->c = f;
        g->node = a3_alloc(a, (size_t)(f.nx + 1) * (f.ny + 1) * (f.ns + 1), sizeof(V3));
        g->cell = a3_alloc(a, (size_t)f.nx * f.ny * f.ns, sizeof(TubeCell));
        for (int axis = 0; axis < 3; axis++)
            g->face[axis] = a3_alloc(a, tube_face_count(g, axis), sizeof(TubeFace));
    }
    double stats[] = {g->max_alignment_error, g->max_parallel_electric_metric_error, g->min_volume,
                      g->max_magnetic_flux_jump};
    void *arrays[] = {stats, g->node, g->cell, g->face[0], g->face[1], g->face[2]};
    size_t bytes[] = {sizeof(stats),
                      (size_t)(f.nx + 1) * (f.ny + 1) * (f.ns + 1) * sizeof(V3),
                      (size_t)f.nx * f.ny * f.ns * sizeof(TubeCell),
                      (size_t)tube_face_count(g, 0) * sizeof(TubeFace),
                      (size_t)tube_face_count(g, 1) * sizeof(TubeFace),
                      (size_t)tube_face_count(g, 2) * sizeof(TubeFace)};
    uint64_t checksum = UINT64_C(14695981039346656037);
    for (int v = 0; v < 6; v++) {
        if ((write ? fwrite(arrays[v], 1, bytes[v], stream)
                   : fread(arrays[v], 1, bytes[v], stream)) != bytes[v])
            a3_fail(a, "geometry cache payload size");
        checksum = hash(checksum, arrays[v], bytes[v]);
    }
    if (write) {
        if (fwrite(&checksum, 1, sizeof(checksum), stream) != sizeof(checksum) || fclose(stream) ||
            rename(temporary, c.geometry_cache))
            a3_fail(a, "geometry cache commit");
    } else {
        uint64_t expected;
        if (fread(&expected, 1, sizeof(expected), stream) != sizeof(expected) ||
            checksum != expected || fgetc(stream) != EOF || ferror(stream))
            a3_fail(a, "geometry cache payload checksum");
        fclose(stream);
        g->max_alignment_error = stats[0];
        g->max_parallel_electric_metric_error = stats[1];
        g->min_volume = stats[2];
        g->max_magnetic_flux_jump = stats[3];
    }
    fprintf(stdout, "Geometry cache %s: %s\n", write ? "created" : "verified and loaded",
            c.geometry_cache);
    fflush(stdout);
    return 1;
}
