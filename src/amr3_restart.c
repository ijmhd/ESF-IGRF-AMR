/*
 * Legacy ESFAMR01 portable checkpoint and dispatch to structure-aware versions. Integer and
 * floating data are encoded little-endian and bound to input signatures.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Portable little-endian field-aligned AMR checkpoint. Ownership is reconstructed
 * from the saved mask, independently of rank count. Quiet counters and travel/
 * regrid clocks preserve hysteresis and the next refinement decision.
 */
#include "amr3.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    A3 *a;
    FILE *f;
    uint64_t hash;
    int write;
} Stream;
/**
 * Read or write exactly the requested bytes while updating the checkpoint checksum; short
 * transfers are fatal.
 */
static void bytes(Stream *s, void *p, size_t n) {
    if ((s->write ? fwrite(p, 1, n, s->f) : fread(p, 1, n, s->f)) != n)
        a3_fail(s->a, "AMR checkpoint transfer size");
    unsigned char *b = p;
    for (size_t k = 0; k < n; k++) {
        s->hash ^= b[k];
        s->hash *= UINT64_C(1099511628211);
    }
}
/**
 * Encode/decode an unsigned 64-bit integer in little-endian order and update the stream
 * checksum.
 */
static uint64_t integer(Stream *s, uint64_t value) {
    unsigned char b[8];
    for (int k = 0; k < 8; k++)
        b[k] = (unsigned char)(value >> (8 * k));
    bytes(s, b, 8);
    value = 0;
    for (int k = 0; k < 8; k++)
        value |= (uint64_t)b[k] << (8 * k);
    return value;
}
/**
 * Transfer floating data in portable little-endian order, retaining the checksum over the actual
 * encoded bytes.
 */
static void reals(Stream *s, double *p, size_t n) {
    unsigned short one = 1;
    if (*(unsigned char *)&one)
        bytes(s, p, n * sizeof(double));
    else
        for (size_t k = 0; k < n; k++) {
            uint64_t v;
            memcpy(&v, p + k, 8);
            v = integer(s, v);
            memcpy(p + k, &v, 8);
        }
}
/**
 * Build the ordered physical/numerical restart signature; changing a signed parameter
 * intentionally rejects incompatible restart.
 */
static void signature(Config3 c, double values[41]) {
    if (c.seed_offset_deg != 0) {
        fprintf(stderr,
                "Legacy AMR checkpoints omit seed longitude; use structure-AMR checkpoints\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    solver3_checkpoint_parameters(&c, values);
    double tail[] = {c.amr_mode,       c.amr_case,        c.amr_interval, c.amr_buffer,  c.amr_hold,
                     c.amr_deplete_on, c.amr_deplete_off, c.amr_grad_on,  c.amr_grad_off};
    memcpy(values + 32, tail, sizeof(tail));
}
/**
 * Write or verify the final checksum, reject trailing data on reads, and close the checkpoint
 * stream.
 */
static void finish(Stream *s) {
    uint64_t hash = s->hash;
    if (integer(s, hash) != hash)
        a3_fail(s->a, "AMR checkpoint checksum");
    if (!s->write && (fgetc(s->f) != EOF || ferror(s->f)))
        a3_fail(s->a, "AMR checkpoint trailing data");
    if (fclose(s->f))
        a3_fail(s->a, "AMR checkpoint close");
}
/**
 * Write ESFAMR01 for legacy masks or dispatch to the structure-aware writer; retain little-
 * endian encoding and atomic replacement.
 */
void a3_checkpoint_write(A3 *a, const double *n, const double *m, const double *marker, double time,
                         int step, double mass0, const int *quiet, double next_regrid,
                         double travel, int retries) {
    if (a->c.amr_structure) {
        a3_checkpoint2_write(a, n, m, marker, time, step, mass0, next_regrid, travel, retries);
        return;
    }
    Solver3 status = {.rank = a->rank};
    uint64_t header[] = {a->nx,
                         a->ny,
                         a->ns,
                         a->nleaf,
                         step,
                         a->solves,
                         retries,
                         a->mesh_epoch,
                         a->regrids,
                         a->refinements,
                         a->coarsenings,
                         41,
                         solver3_file_hash(&status, a->c.coeff),
                         solver3_file_hash(&status, a->c.collision_profile)};
    double info[] = {time,
                     mass0,
                     a->boundary_total,
                     a->max_residual,
                     a->remap_error,
                     next_regrid,
                     travel,
                     a->import_mass_error,
                     a->geometry_density_adjustment,
                     a->max_potential_remap_change,
                     a->material_min,
                     a->material_max};
    char path[1200], part[1200];
    snprintf(path, sizeof(path), "%s/checkpoint_latest.bin", a->c.out);
    snprintf(part, sizeof(part), "%s/checkpoint_latest.part", a->c.out);
    Stream s = {a, fopen(part, "wb"), UINT64_C(14695981039346656037), 1};
    if (!s.f)
        a3_fail(a, "AMR checkpoint open");
    char magic[] = "ESFAMR01";
    bytes(&s, magic, 8);
    for (int k = 0; k < 14; k++)
        integer(&s, header[k]);
    reals(&s, info, 12);
    double params[41];
    signature(a->c, params);
    reals(&s, params, 41);
    char code[] = ESF_GEOMETRY_ID;
    bytes(&s, code, 128);
    bytes(&s, a->refined, a->nb);
    for (int b = 0; b < a->nb; b++)
        integer(&s, quiet[b]);
    reals(&s, (double *)n, (size_t)a->nleaf * a->ns);
    reals(&s, (double *)m, (size_t)a->nleaf * a->ns);
    size_t nv = (size_t)(a->fx + 1) * (a->fy + 1);
    double *phi = a->phi ? a->phi : a3_alloc(a, nv, sizeof(double));
    reals(&s, phi, nv);
    if (!a->phi)
        free(phi);
    finish(&s);
    if (rename(part, path))
        a3_fail(a, "AMR checkpoint atomic rename");
}
/**
 * Recognize the checkpoint magic, validate legacy signatures or dispatch to the structure
 * reader, then scatter conserved states to current rank ownership.
 */
int a3_checkpoint_read(A3 *a, TubeGrid *g, Config3 c, int rank, int size, double **n, double **m,
                       double *time, int *step, double *mass0, int **quiet, double *next_regrid,
                       double *travel, int *retries) {
    A3 status = {.rank = rank, .size = size};
    int native = 0;
    if (!rank && c.restart[0]) {
        FILE *f = fopen(c.restart, "rb");
        if (!f)
            a3_fail(&status, "restart open");
        char magic[8];
        if (fread(magic, 1, 8, f) == 8) {
            if (!memcmp(magic, "ESFAMR01", 8))
                native = 1;
            else if (!memcmp(magic, "ESFAMR02", 8) || !memcmp(magic, "ESFAMR03", 8) ||
                     (!memcmp(magic, "ESFAMR04", 8) || !memcmp(magic, "ESFAMR05", 8)))
                native = 2;
        }
        fclose(f);
    }
    MPI_Bcast(&native, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!native)
        return 0;
    if (native == 2)
        return a3_checkpoint2_read(a, g, c, rank, size, n, m, time, step, mass0, quiet, next_regrid,
                                   travel, retries);
    if (c.amr_structure || c.amr_slab)
        a3_fail(
            &status,
            "legacy AMR checkpoints do not contain local masks, physical-time hysteresis or protected-mode markers");
    uint64_t header[14] = {0};
    double info[12] = {0};
    int nb = c.grid.nx * c.grid.ny;
    unsigned char *mask = a3_alloc(&status, nb, 1);
    *quiet = a3_alloc(&status, nb, sizeof(int));
    double *N = NULL, *M = NULL,
           *phi =
               a3_alloc(&status, (size_t)(2 * c.grid.nx + 1) * (2 * c.grid.ny + 1), sizeof(double));
    if (!rank) {
        Solver3 check = {.rank = rank};
        Stream s = {&status, fopen(c.restart, "rb"), UINT64_C(14695981039346656037), 0};
        if (!s.f)
            a3_fail(&status, "AMR restart open");
        char magic[8];
        bytes(&s, magic, 8);
        for (int k = 0; k < 14; k++)
            header[k] = integer(&s, 0);
        reals(&s, info, 12);
        if (header[0] != (uint64_t)c.grid.nx || header[1] != (uint64_t)c.grid.ny ||
            header[2] != (uint64_t)c.grid.ns || header[3] < 1 || header[3] > 4 * (uint64_t)nb ||
            header[11] != 41)
            a3_fail(&status, "AMR restart grid/version mismatch");
        for (int k = 4; k <= 10; k++)
            if (header[k] > INT_MAX)
                a3_fail(&status, "AMR restart counter overflow");
        if (header[12] != solver3_file_hash(&check, c.coeff) ||
            header[13] != solver3_file_hash(&check, c.collision_profile))
            a3_fail(&status, "AMR restart input hash mismatch");
        for (int k = 0; k < 10; k++)
            if (!isfinite(info[k]))
                a3_fail(&status, "AMR restart nonfinite history");
        if (info[0] < 0 || info[0] > c.end || info[1] <= 0 || info[3] < 0 ||
            (c.amr_mode == 3 && info[5] < info[0] - 1e-9) || info[6] < 0 || isnan(info[10]) ||
            isnan(info[11]))
            a3_fail(&status, "AMR restart history mismatch");
        double params[41], expected[41];
        reals(&s, params, 41);
        signature(c, expected);
        for (int k = 0; k < 41; k++)
            if (!isfinite(params[k]) ||
                fabs(params[k] - expected[k]) > 5e-14 * fmax(1., fabs(expected[k])))
                a3_fail(&status, "AMR restart physical/numerical/refinement mismatch");
        char code[128];
        bytes(&s, code, 128);
        if (!a3_geometry_identity(&status, code))
            a3_fail(&status, "AMR restart geometry implementation mismatch");
        bytes(&s, mask, nb);
        for (int b = 0; b < nb; b++) {
            uint64_t q = integer(&s, 0);
            if (mask[b] > 1 || q > (uint64_t)c.amr_hold)
                a3_fail(&status, "AMR restart mask/quiet mismatch");
            (*quiet)[b] = (int)q;
        }
        size_t cells = header[3] * c.grid.ns;
        N = a3_alloc(&status, cells, sizeof(double));
        M = a3_alloc(&status, cells, sizeof(double));
        reals(&s, N, cells);
        reals(&s, M, cells);
        reals(&s, phi, (size_t)(2 * c.grid.nx + 1) * (2 * c.grid.ny + 1));
        finish(&s);
        for (size_t k = 0; k < cells; k++)
            if (!isfinite(N[k]) || !isfinite(M[k]) || N[k] < 0 || (N[k] == 0 && M[k] != 0))
                a3_fail(&status, "AMR restart invalid state");
        for (int k = 0; k < (2 * c.grid.nx + 1) * (2 * c.grid.ny + 1); k++)
            if (!isfinite(phi[k]))
                a3_fail(&status, "AMR restart invalid potential");
    }
    MPI_Bcast(header, 14, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(info, 12, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(mask, nb, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(*quiet, nb, MPI_INT, 0, MPI_COMM_WORLD);
    a3_mesh_mask(a, g, c, rank, size, c.amr_mode, mask, NULL, NULL);
    free(mask);
    if (a->nleaf != (int)header[3])
        a3_fail(a, "AMR restart leaf count mismatch");
    *n = a3_alloc(a, (size_t)a->slots * a->ns, sizeof(double));
    *m = a3_alloc(a, (size_t)a->slots * a->ns, sizeof(double));
    int *counts = a3_alloc(a, 2 * size, sizeof(int)), *offset = counts + size;
    for (int r = 0; r < size; r++) {
        offset[r] = a->row[a->partition[r]] * a->ns;
        counts[r] = (a->row[a->partition[r + 1]] - a->row[a->partition[r]]) * a->ns;
    }
    MPI_Scatterv(N, counts, offset, MPI_DOUBLE, *n, a->owned_count * a->ns, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    MPI_Scatterv(M, counts, offset, MPI_DOUBLE, *m, a->owned_count * a->ns, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    free(N);
    free(M);
    free(counts);
    MPI_Bcast(phi, (a->fx + 1) * (a->fy + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (c.amr_case != 1) {
        a3_potential_create(a);
        a3_potential_seed(a, phi);
    }
    free(phi);
    *time = info[0];
    *mass0 = info[1];
    a->boundary_total = info[2];
    a->max_residual = info[3];
    a->remap_error = info[4];
    *next_regrid = info[5];
    *travel = info[6];
    a->import_mass_error = info[7];
    a->geometry_density_adjustment = info[8];
    a->max_potential_remap_change = info[9];
    a->material_min = info[10];
    a->material_max = info[11];
    *step = (int)header[4];
    a->solves = (int)header[5];
    *retries = (int)header[6];
    a->mesh_epoch = (int)header[7];
    a->regrids = (int)header[8];
    a->refinements = (int)header[9];
    a->coarsenings = (int)header[10];
    if (!rank) {
        fprintf(stdout, "AMR restart t=%.17g step=%d leaves=%d on %d MPI ranks\n", *time, *step,
                a->nleaf, size);
        fflush(stdout);
    }
    return 1;
}
