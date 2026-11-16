/*
 * Structure-aware ESFAMR02/03/04/05 checkpoint formats. Preserve density/field masks, hysteresis
 * clocks, optional seed tracer and conservation history.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Structure-AMR restart includes density masks, physical clocks and tracer.
 */
#include "amr3.h"
#include <math.h>
#include <limits.h>
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
        a3_fail(s->a, "structure checkpoint size");
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
 * Write or verify the final checksum, reject trailing data on reads, and close the checkpoint
 * stream.
 */
static void finish(Stream *s) {
    uint64_t hash = s->hash;
    if (integer(s, hash) != hash)
        a3_fail(s->a, "structure checkpoint checksum");
    if (!s->write && (fgetc(s->f) != EOF || ferror(s->f)))
        a3_fail(s->a, "structure checkpoint trailing bytes");
    if (fclose(s->f))
        a3_fail(s->a, "structure checkpoint close");
}
/**
 * Choose the exact restart signature length required by this case's driver/neutral extensions.
 */
static int parameters(Config3 c) {
    return c.spatial_neutral ? 72 : (c.driver || c.background_east != 0 ? 68 : 64);
}
/**
 * Build the ordered physical/numerical restart signature; changing a signed parameter
 * intentionally rejects incompatible restart.
 */
static void signature(Config3 c, double value[72]) {
    memset(value, 0, 72 * sizeof(double));
    solver3_checkpoint_parameters(&c, value);
    double tail[] = {
        c.amr_mode,           c.amr_case,        c.amr_interval,    c.amr_buffer,
        c.amr_hold,           c.amr_deplete_on,  c.amr_deplete_off, c.amr_grad_on,
        c.amr_grad_off,       c.amr_structure,   c.amr_slab,        c.amr_seed_guard,
        c.amr_fixture,        c.amr_curve_on,    c.amr_curve_off,   c.amr_velocity_on,
        c.amr_sigma_on,       c.amr_quiet_s,     c.amr_seed_lambda, c.amr_branch_lambda,
        c.amr_seed_amplitude, c.amr_work_weight, c.amr_coo,         c.amr_residual};
    memcpy(value + 32, tail, sizeof(tail));
    value[56] = c.seed_offset_deg; /* Reserved signature slot; legacy centered seeds stored zero. */
    if (c.night_chemistry) {
        value[57] = 1;
        value[58] = c.utc_start;
        value[59] = c.chem_Ti;
        value[60] = c.chem_Te;
        value[61] = c.night_reactions;
        value[62] = 2;
    }
    if (parameters(c) >= 68) {
        uint64_t h = c.driver ? c.driver->hash : 0;
        value[63] = c.driver || c.background_east != 0;
        value[64] = (double)(h >> 32);
        value[65] = (double)(h & UINT64_C(0xffffffff));
        value[66] = 1;
    }
    if (c.amr_levels == 3)
        value[67] = 3;
    if (c.spatial_neutral) {
        value[68] = 1;
        value[69] = (double)(c.neutral_grid_hash >> 32);
        value[70] = (double)(c.neutral_grid_hash & UINT64_C(0xffffffff));
        value[71] = 1;
    }
}
/**
 * Atomically write structure-aware masks/history, optional tracer, actual conserved states,
 * potential and checksummed physical/numerical signatures.
 */
void a3_checkpoint2_write(A3 *a, const double *n, const double *m, const double *marker,
                          double time, int step, double mass0, double next_regrid, double travel,
                          int retries) {
    Solver3 check = {.rank = a->rank};
    uint64_t h[] = {a->nx,
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
                    parameters(a->c),
                    solver3_file_hash(&check, a->c.coeff),
                    solver3_file_hash(&check, a->c.collision_profile)};
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
        a3_fail(a, "structure checkpoint open");
    char magic[9];
    strcpy(magic, a->c.amr_levels == 3     ? "ESFAMR05"
                  : parameters(a->c) == 72 ? "ESFAMR04"
                                           : (parameters(a->c) == 68 ? "ESFAMR03" : "ESFAMR02"));
    char envelope[] = "ESFSRC01", source[] = ESF_NIGHT_ID;
    bytes(&s, envelope, 8);
    bytes(&s, source, 64);
    bytes(&s, magic, 8);
    for (int k = 0; k < 14; k++)
        integer(&s, h[k]);
    reals(&s, info, 12);
    double params[72];
    signature(a->c, params);
    reals(&s, params, parameters(a->c));
    char code[] = ESF_GEOMETRY_ID;
    bytes(&s, code, 128);
    bytes(&s, a->refined, a->nb);
    bytes(&s, a->density_refined, (size_t)a->nb * a->ns);
    bytes(&s, a->field_requested, a->nb);
    reals(&s, a->quiet_since, a->mask_count);
    reals(&s, a->born, a->mask_count);
    integer(&s, marker != NULL);
    if (a->c.amr_levels == 3) {
        size_t nc = (size_t)a->nb / 4 * a->ns;
        bytes(&s, a->macro_mask, nc);
        reals(&s, a->macro_quiet, nc);
        reals(&s, a->macro_born, nc);
    }
    size_t cells = (size_t)a->nleaf * a->ns;
    reals(&s, (double *)n, cells);
    reals(&s, (double *)m, cells);
    if (marker)
        reals(&s, (double *)marker, cells);
    size_t nv = (size_t)(a->fx + 1) * (a->fy + 1);
    double *phi = a->phi ? a->phi : a3_alloc(a, nv, sizeof(double));
    reals(&s, phi, nv);
    if (!a->phi)
        free(phi);
    finish(&s);
    if (rename(part, path))
        a3_fail(a, "structure checkpoint rename");
}
/**
 * Validate and reconstruct a structure-aware checkpoint, distribute states under the current MPI
 * decomposition, and restore mesh/clocks/history.
 */
int a3_checkpoint2_read(A3 *a, TubeGrid *g, Config3 c, int rank, int size, double **n, double **m,
                        double *time, int *step, double *mass0, int **quiet, double *next_regrid,
                        double *travel, int *retries) {
    A3 status = {.rank = rank, .size = size};
    int nb = c.grid.nx * c.grid.ny, nmask = nb * (c.amr_slab ? c.grid.ns : 1);
    if (!c.amr_structure)
        a3_fail(&status, "structure checkpoint requires structure criterion");
    uint64_t h[14] = {0}, has_marker = 0;
    double info[12] = {0};
    unsigned char *field = a3_alloc(&status, nb, 1),
                  *density = a3_alloc(&status, (size_t)nb * c.grid.ns, 1),
                  *requested = a3_alloc(&status, nb, 1);
    unsigned char *macro = a3_alloc(&status, (size_t)nb / 4 * c.grid.ns, 1);
    double *mq = a3_alloc(&status, (size_t)nb / 4 * c.grid.ns, sizeof(double)),
           *mb = a3_alloc(&status, (size_t)nb / 4 * c.grid.ns, sizeof(double));
    double *since = a3_alloc(&status, nmask, sizeof(double)),
           *born = a3_alloc(&status, nmask, sizeof(double));
    size_t nv = (size_t)(2 * c.grid.nx + 1) * (2 * c.grid.ny + 1);
    double *phi = a3_alloc(&status, nv, sizeof(double)), *N = NULL, *M = NULL, *P = NULL;
    if (!rank) {
        Solver3 check = {.rank = rank};
        Stream s = {&status, fopen(c.restart, "rb"), UINT64_C(14695981039346656037), 0};
        if (!s.f)
            a3_fail(&status, "structure restart open");
        char envelope[8], source[64], magic[8];
        bytes(&s, envelope, 8);
        bytes(&s, source, 64);
        if (memcmp(envelope, "ESFSRC01", 8) || memcmp(source, ESF_NIGHT_ID, 64))
            a3_fail(&status, "AMR restart numerical source identity mismatch");
        bytes(&s, magic, 8);
        for (int k = 0; k < 14; k++)
            h[k] = integer(&s, 0);
        reals(&s, info, 12);
        if (memcmp(magic,
                   c.amr_levels == 3     ? "ESFAMR05"
                   : parameters(c) == 72 ? "ESFAMR04"
                                         : (parameters(c) == 68 ? "ESFAMR03" : "ESFAMR02"),
                   8) ||
            h[0] != (uint64_t)c.grid.nx || h[1] != (uint64_t)c.grid.ny ||
            h[2] != (uint64_t)c.grid.ns || h[3] < (uint64_t)nb || h[3] > 4 * (uint64_t)nb ||
            h[11] != (uint64_t)parameters(c))
            a3_fail(&status, "structure restart layout/version");
        for (int k = 4; k <= 10; k++)
            if (h[k] > INT_MAX)
                a3_fail(&status, "structure restart counter range");
        if (h[12] != solver3_file_hash(&check, c.coeff) ||
            h[13] != solver3_file_hash(&check, c.collision_profile))
            a3_fail(&status, "structure restart input hashes");
        for (int k = 0; k < 10; k++)
            if (!isfinite(info[k]))
                a3_fail(&status, "structure restart nonfinite history");
        if (isnan(info[10]) || isnan(info[11]))
            a3_fail(&status, "structure restart invalid material bounds");
        if (info[0] < 0 || info[0] > c.end || info[1] <= 0 ||
            (info[5] < 0 || (c.amr_mode == 3 && info[5] < info[0] - 1e-9)) || info[6] < 0)
            a3_fail(&status, "structure restart clocks");
        double params[72], expected[72];
        reals(&s, params, parameters(c));
        signature(c, expected);
        for (int k = 0; k < parameters(c); k++)
            if (!isfinite(params[k]) ||
                fabs(params[k] - expected[k]) > 5e-14 * fmax(1., fabs(expected[k])))
                a3_fail(&status, "structure restart parameters changed");
        char code[128];
        bytes(&s, code, 128);
        if (!a3_geometry_identity(&status, code))
            a3_fail(&status, "structure restart geometry identity");
        bytes(&s, field, nb);
        bytes(&s, density, (size_t)nb * c.grid.ns);
        bytes(&s, requested, nb);
        reals(&s, since, nmask);
        reals(&s, born, nmask);
        has_marker = integer(&s, 0);
        if (has_marker > 1)
            a3_fail(&status, "structure restart tracer flag");
        if (c.amr_levels == 3) {
            size_t nc = (size_t)nb / 4 * c.grid.ns;
            bytes(&s, macro, nc);
            reals(&s, mq, nc);
            reals(&s, mb, nc);
            for (size_t z = 0; z < nc; z++)
                if (macro[z] > 1 || !isfinite(mq[z]) || !isfinite(mb[z]))
                    a3_fail(&status, "invalid macro checkpoint");
        }
        for (int b = 0; b < nb; b++) {
            if (field[b] > 1 || requested[b] > 1)
                a3_fail(&status, "invalid field mask");
            for (int k = 0; k < c.grid.ns; k++)
                if (density[b * c.grid.ns + k] > 1 || (density[b * c.grid.ns + k] && !field[b]))
                    a3_fail(&status, "invalid local density mask");
        }
        for (int q = 0; q < nmask; q++)
            if (!isfinite(since[q]) || !isfinite(born[q]) || since[q] < -1 || since[q] > info[0] ||
                born[q] < 0 || born[q] > info[0])
                a3_fail(&status, "invalid physical hysteresis history");
        size_t cells = h[3] * c.grid.ns;
        N = a3_alloc(&status, cells, sizeof(double));
        M = a3_alloc(&status, cells, sizeof(double));
        if (has_marker)
            P = a3_alloc(&status, cells, sizeof(double));
        reals(&s, N, cells);
        reals(&s, M, cells);
        if (P)
            reals(&s, P, cells);
        reals(&s, phi, nv);
        finish(&s);
        for (size_t q = 0; q < cells; q++)
            if (!isfinite(N[q]) || N[q] < 0 || !isfinite(M[q]) || (P && !isfinite(P[q])))
                a3_fail(&status, "invalid structure state");
    }
    MPI_Bcast(h, 14, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(info, 12, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&has_marker, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(field, nb, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(density, nb * c.grid.ns, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(requested, nb, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(since, nmask, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(born, nmask, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    unsigned char *backing = a3_alloc(&status, nmask, 1);
    for (int b = 0; b < nb; b++)
        for (int k = 0; k < (c.amr_slab ? c.grid.ns : 1); k++)
            backing[b * (c.amr_slab ? c.grid.ns : 1) + k] = field[b];
    a3_mesh_mask_work(a, g, c, rank, size, c.amr_mode, backing, NULL, NULL, density);
    free(backing);
    if (a->nleaf != (int)h[3])
        a3_fail(a, "structure restart backing count");
    if (c.amr_levels == 3) {
        int nc = nb / 4 * c.grid.ns;
        MPI_Bcast(macro, nc, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(mq, nc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(mb, nc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        memcpy(a->macro_mask, macro, nc);
        memcpy(a->macro_next, macro, nc);
        memcpy(a->macro_quiet, mq, nc * sizeof(double));
        memcpy(a->macro_born, mb, nc * sizeof(double));
    }
    free(macro);
    free(mq);
    free(mb);
    memcpy(a->density_refined, density, (size_t)nb * c.grid.ns);
    memcpy(a->field_requested, requested, nb);
    a3_segment_metrics(a);
    a3_face_topology(a);
    memcpy(a->quiet_since, since, nmask * sizeof(double));
    memcpy(a->born, born, nmask * sizeof(double));
    free(field);
    free(density);
    free(requested);
    free(since);
    free(born);
    *n = a3_alloc(a, (size_t)a->slots * a->ns, sizeof(double));
    *m = a3_alloc(a, (size_t)a->slots * a->ns, sizeof(double));
    if (has_marker)
        a->marker = a3_alloc(a, (size_t)a->slots * a->ns, sizeof(double));
    int *counts = a3_alloc(a, 2 * size, sizeof(int)), *offset = counts + size;
    for (int r = 0; r < size; r++) {
        offset[r] = a->row[a->partition[r]] * a->ns;
        counts[r] = (a->row[a->partition[r + 1]] - a->row[a->partition[r]]) * a->ns;
    }
    MPI_Scatterv(N, counts, offset, MPI_DOUBLE, *n, a->owned_count * a->ns, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    MPI_Scatterv(M, counts, offset, MPI_DOUBLE, *m, a->owned_count * a->ns, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    if (has_marker)
        MPI_Scatterv(P, counts, offset, MPI_DOUBLE, a->marker, a->owned_count * a->ns, MPI_DOUBLE,
                     0, MPI_COMM_WORLD);
    free(N);
    free(M);
    free(P);
    free(counts);
    MPI_Bcast(phi, (int)nv, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (c.amr_case != 1) {
        a3_potential_create(a);
        a3_potential_seed(a, phi);
    }
    free(phi);
    *quiet = a3_alloc(a, nmask, sizeof(int));
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
    *step = (int)h[4];
    a->solves = (int)h[5];
    *retries = (int)h[6];
    a->mesh_epoch = (int)h[7];
    a->regrids = (int)h[8];
    a->refinements = (int)h[9];
    a->coarsenings = (int)h[10];
    a->tag_time = *time;
    if (!a3_valid(a, *n, *m))
        a3_fail(a, "invalid protected marker or density after restart");
    return 1;
}
