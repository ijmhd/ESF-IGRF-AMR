/*
 * Uniform-grid portable checkpoint reader and physical/numerical input signature. It is distinct
 * from native structure-AMR restart.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* Portable little-endian accepted-state restart; ETD2 is a one-step method. */
#include "solver.h"
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/**
 * Read one little-endian 64-bit unsigned checkpoint word; short reads are fatal.
 */
static uint64_t read64(Solver3 *s, FILE *f) {
    unsigned char bytes[8];
    if (fread(bytes, 1, 8, f) != 8)
        fail3(s, "truncated restart header");
    uint64_t value = 0;
    for (int i = 0; i < 8; i++)
        value |= (uint64_t)bytes[i] << (8 * i);
    return value;
}
/**
 * Decode one checkpoint IEEE double from its little-endian bit pattern.
 */
static double real64(Solver3 *s, FILE *f) {
    uint64_t bits = read64(s, f);
    double x;
    memcpy(&x, &bits, 8);
    return x;
}
/**
 * Calculate an exact input-file FNV identity, with zero for an absent optional input.
 */
uint64_t solver3_file_hash(Solver3 *s, const char *path) {
    if (!path[0])
        return 0;
    FILE *f = fopen(path, "rb");
    if (!f)
        fail3(s, "cannot open restart input for identity check");
    uint64_t h = UINT64_C(14695981039346656037);
    unsigned char block[8192];
    size_t n;
    while ((n = fread(block, 1, sizeof(block), f)))
        for (size_t i = 0; i < n; i++) {
            h ^= block[i];
            h *= UINT64_C(1099511628211);
        }
    if (ferror(f))
        fail3(s, "cannot hash restart input");
    fclose(f);
    return h;
}
/**
 * Read a portable array of doubles into caller storage.
 */
static void array64(Solver3 *s, FILE *f, double *x, int n) {
    unsigned short one = 1;
    if (*(unsigned char *)&one) {
        if (fread(x, sizeof(double), n, f) != (size_t)n)
            fail3(s, "truncated restart state");
    } else
        for (int i = 0; i < n; i++)
            x[i] = real64(s, f);
    for (int i = 0; i < n; i++)
        if (!isfinite(x[i]))
            fail3(s, "nonfinite restart state");
}
/**
 * Construct the fixed-order uniform-grid physical/numerical signature for compatibility checks.
 */
void solver3_checkpoint_parameters(const Config3 *c, double expected[32]) {
    double values[] = {c->grid.lon,
                       c->grid.span,
                       c->grid.hmin,
                       c->grid.hmax,
                       c->grid.hbase,
                       c->grid.trace_step,
                       c->degree,
                       c->aligned_dipole,
                       c->year,
                       c->dipole_equator_nT,
                       c->parallel,
                       c->hall,
                       c->temperature,
                       c->periodic_x,
                       c->open_apex,
                       c->open_ends,
                       c->parallel_order,
                       c->parallel_relaxation,
                       c->exb_pdm4,
                       c->pdm_beta,
                       c->background_profile,
                       c->background_east,
                       c->regional_x,
                       c->cfl,
                       c->dtmax,
                       c->rtol,
                       c->amplitude,
                       c->seed_height,
                       c->seed_width,
                       c->seed_east,
                       c->e_peak,
                       2};
    memcpy(expected, values, sizeof(values));
}
/**
 * Restore and validate the uniform-grid checkpoint, including original conservation history and
 * imposed-field compatibility.
 */
void solver3_restart(Solver3 *s, const char *path, double *time, int *step) {
    if (s->c.seed_offset_deg != 0)
        fail3(
            s,
            "legacy static restart has no seed longitude identity; use native structure-AMR restart");
    double info[4] = {0};
    int counts[3] = {0};
    const Config3 *c = &s->c;
    if (!s->rank) {
        FILE *f = fopen(path, "rb");
        if (!f)
            fail3(s, "cannot open restart");
        char magic[8];
        if (fread(magic, 1, 8, f) != 8)
            fail3(s, "truncated restart magic");
        int version = !memcmp(magic, "ESFRST02", 8) ? 2 : (!memcmp(magic, "ESFRST01", 8) ? 1 : 0);
        if (!version)
            fail3(s, "invalid restart version");
        if (version == 1 && (!c->periodic_x || c->e_peak != 0))
            fail3(s, "V1 restart omits E-peak and regional boundary model; require ESFRST02");
        uint64_t nx = read64(s, f), ny = read64(s, f), ns = read64(s, f);
        if (nx != (uint64_t)c->grid.nx || ny != (uint64_t)c->grid.ny || ns != (uint64_t)c->grid.ns)
            fail3(s, "restart grid mismatch");
        for (int i = 0; i < 3; i++) {
            uint64_t value = read64(s, f);
            if (value > INT_MAX)
                fail3(s, "restart counter overflow");
            counts[i] = (int)value;
        }
        int parameters = version == 2 ? 32 : 30;
        if (read64(s, f) != (uint64_t)parameters)
            fail3(s, "restart parameter version mismatch");
        if (read64(s, f) != solver3_file_hash(s, c->coeff) ||
            read64(s, f) != solver3_file_hash(s, c->collision_profile))
            fail3(s, "restart background/field input mismatch");
        for (int i = 0; i < 4; i++)
            info[i] = real64(s, f);
        if (!isfinite(info[0]) || info[0] < 0 || info[0] > c->end || !isfinite(info[1]) ||
            info[1] <= 0 || !isfinite(info[2]) || !isfinite(info[3]) || info[3] < 0)
            fail3(s, "invalid restart time or conservation history");
        double expected[32];
        solver3_checkpoint_parameters(c, expected);
        for (int i = 0; i < parameters; i++) {
            double value = real64(s, f);
            if (!isfinite(value) ||
                fabs(value - expected[i]) > 5e-14 * fmax(1., fabs(expected[i]))) {
                fprintf(stderr, "restart parameter %d: stored %.17g requested %.17g\n", i, value,
                        expected[i]);
                fail3(s, "restart physical/numerical configuration mismatch");
            }
        }
        array64(s, f, s->rho, s->total);
        array64(s, f, s->momentum, s->total);
        array64(s, f, s->phi, (c->grid.nx + 1) * (c->grid.ny + 1));
        if (fgetc(f) != EOF || ferror(f))
            fail3(s, "restart file size mismatch");
        fclose(f);
    }
#ifdef ESF_MPI
    MPI_Bcast(info, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(counts, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(s->rho, s->counts, s->offsets, MPI_DOUBLE,
                 s->rank ? s->rho + s->first * c->grid.ns : MPI_IN_PLACE, s->n, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    MPI_Scatterv(s->momentum, s->counts, s->offsets, MPI_DOUBLE,
                 s->rank ? s->momentum + s->first * c->grid.ns : MPI_IN_PLACE, s->n, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    MPI_Bcast(s->phi, (c->grid.nx + 1) * (c->grid.ny + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    *time = info[0];
    s->mass0 = info[1];
    s->boundary_total = info[2];
    s->max_potential_residual = info[3];
    *step = counts[0];
    s->solves = counts[1];
    s->retries = counts[2];
    if (!solver3_valid(s, s->rho, s->momentum))
        fail3(s, "invalid restarted density/momentum");
    if (s->p.phi)
        for (int row = 0; row < s->p.n; row++)
            s->p.phi[row] = s->phi[row] - s->background_potential_gradient *
                                              ((double)(row % (c->grid.nx + 1)) / c->grid.nx - .5);
#ifdef ESF_PETSC
    if (c->petsc_potential && solver3_petsc_seed(s))
        fail3(s, "cannot seed restarted potential");
#endif
    if (!s->rank)
        fprintf(stdout,
                "Restart loaded at t=%.17g, step=%d; original mass and boundary history retained\n",
                *time, *step);
}
