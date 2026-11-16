/*
 * Collective mesh rebuild, ownership migration and conservative transfer. Geometry remains
 * fixed; density and electric-potential approximation spaces may change separately.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only */
#include "amr3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * Transfer O+ density/momentum/optional marker intersections between old and new rank-owned
 * slabs before conservative remapping.
 */
static void migrate(A3 *old, A3 *moved, double *n, double *m, double *N, double *M) {
    int *s = a3_alloc(old, 4 * old->size, sizeof(int)), *so = s + old->size, *r = so + old->size,
        *ro = r + old->size;
    for (int peer = 0; peer < old->size; peer++) {
        int lo = old->row_first > moved->partition[peer] ? old->row_first : moved->partition[peer];
        int hi =
            old->row_last < moved->partition[peer + 1] ? old->row_last : moved->partition[peer + 1];
        if (hi > lo) {
            s[peer] = (old->row[hi] - old->row[lo]) * old->ns;
            so[peer] = (old->row[lo] - old->owned_first) * old->ns;
        }
        lo = moved->row_first > old->partition[peer] ? moved->row_first : old->partition[peer];
        hi =
            moved->row_last < old->partition[peer + 1] ? moved->row_last : old->partition[peer + 1];
        if (hi > lo) {
            r[peer] = (old->row[hi] - old->row[lo]) * old->ns;
            ro[peer] = (old->row[lo] - moved->owned_first) * old->ns;
        }
    }
    MPI_Alltoallv(n, s, so, MPI_DOUBLE, N, r, ro, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Alltoallv(m, s, so, MPI_DOUBLE, M, r, ro, MPI_DOUBLE, MPI_COMM_WORLD);
    if (old->marker)
        MPI_Alltoallv(old->marker, s, so, MPI_DOUBLE, moved->marker, r, ro, MPI_DOUBLE,
                      MPI_COMM_WORLD);
    free(s);
}
/**
 * Detect actual mask/field-space changes, repartition molecules and O+, build/remap the new
 * mesh, verify global invariants and resolve potential. A no-change request returns without
 * rebuilding.
 */
int a3_adapt(A3 *a, double **n, double **m, const unsigned char *mask, double time, int step,
             int exact_initial) {
    int refine = 0, coarsen = 0;
    const unsigned char *oldmask = a3_density_mask(a);
    for (int b = 0; b < a->mask_count; b++) {
        refine += mask[b] && !oldmask[b];
        coarsen += !mask[b] && oldmask[b];
    }
    int field_change = 0;
    for (int b = 0; b < a->nb; b++) {
        int request = a->field_requested[b];
        for (int k = 0; k < a->ns; k++)
            request |= mask[a->c.amr_slab ? b * a->ns + k : b];
        field_change |= request != a->refined[b];
    }
    int macro_change = 0;
    if (a->c.amr_levels == 3)
        macro_change = memcmp(a->macro_mask, a->macro_next, (size_t)a->nb / 4 * a->ns) != 0;
    if (!refine && !coarsen && !field_change && !macro_change)
        return 0;
    /* Plan with the same routine used by mesh construction. Molecular fields
     * retain global backing, but only old/new owned slabs need valid values.
     * Transfer changed ownership before constructors refresh their new halos. */
    int *cuts = a3_alloc(a, a->size + 1, sizeof(int));
    unsigned char *expanded = NULL;
    const unsigned char *work = mask;
    if (!a->c.amr_slab) {
        expanded = a3_alloc(a, (size_t)a->nb * a->ns, 1);
        for (int b = 0; b < a->nb; b++)
            memset(expanded + b * a->ns, mask[b], a->ns);
        work = expanded;
    }
    a3_partition_plan(a, work, cuts);
    free(expanded);
    a3_night_repartition(a, cuts);
    double begin = wall3(), mass = a3_mass(a, *n), momentum = a3_mass(a, *m), norm = 0;
    double marker_mass = a->marker ? a3_mass(a, a->marker) : 0;
    for (int p = 0; p < a->owned_count; p++)
        for (int k = 0; k < a->ns; k++)
            norm += fabs((*m)[p * a->ns + k]) * a->metric[a->global[p] * a->ns + k].volume;
    MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    a3_prepare(a, *n, *m);
    if (a->c.amr_case != 1)
        a3_potential(a);
    A3 next;
    a3_mesh_mask(&next, a->g, a->c, a->rank, a->size, 3, mask, a, cuts);
    free(cuts);
    double *N = a3_alloc(&next, (size_t)next.slots * next.ns, sizeof(double));
    double *M = a3_alloc(&next, (size_t)next.slots * next.ns, sizeof(double));
    if (a->marker)
        next.marker = a3_alloc(&next, (size_t)next.slots * next.ns, sizeof(double));
    const unsigned char *previous = a3_density_mask(a);
    for (int b = 0; b < a->mask_count; b++)
        if (mask[b] && !previous[b])
            next.born[b] = time;
    int migration = memcmp(a->partition, next.partition, (a->size + 1) * sizeof(int)) != 0;
    if (exact_initial) {
        a3_initialize(&next, N, M, 0);
        if (next.marker)
            a3_initialize_marker(&next, N);
        next.boundary_total = a->boundary_total;
    } else if (migration) {
        unsigned char *requested = a->field_requested;
        a->field_requested = a->refined;
        A3 moved;
        a3_mesh_mask(&moved, a->g, a->c, a->rank, a->size, 3, a3_density_mask(a), a,
                     next.partition);
        a->field_requested = requested;
        if (a->c.amr_levels == 3) {
            memcpy(moved.macro_mask, a->macro_mask, (size_t)a->nb / 4 * a->ns);
            a3_segment_metrics(&moved);
            a3_face_topology(&moved);
        }
        if (moved.nleaf != a->nleaf)
            a3_fail(a, "migration must preserve the old backing layout");
        double *mn = a3_alloc(&moved, (size_t)moved.slots * moved.ns, sizeof(double));
        double *mm = a3_alloc(&moved, (size_t)moved.slots * moved.ns, sizeof(double));
        if (a->marker)
            moved.marker = a3_alloc(&moved, (size_t)moved.slots * moved.ns, sizeof(double));
        moved.material_min = a->material_min;
        moved.material_max = a->material_max;
        moved.boundary_total = a->boundary_total;
        migrate(a, &moved, *n, *m, mn, mm);
        a3_remap(&moved, &next, mn, mm, N, M);
        next.t_exchange += moved.t_exchange;
        next.t_prepare += moved.t_prepare;
        free(mn);
        free(mm);
        free(moved.marker);
        a3_free(&moved);
    } else
        a3_remap(a, &next, *n, *m, N, M);
    double particle_error = fabs(a3_mass(&next, N) - mass) / mass;
    double momentum_error = fabs(a3_mass(&next, M) - momentum) / fmax(norm, 1e-300);
    double marker_error =
        next.marker ? fabs(a3_mass(&next, next.marker) - marker_mass) / fmax(marker_mass, 1e-300)
                    : 0;
    if (particle_error > 1e-11 || momentum_error > 1e-11 || marker_error > 1e-11 ||
        !a3_valid(&next, N, M))
        a3_fail(a, "dynamic transfer invariant");
    next.import_mass_error = a->import_mass_error;
    next.geometry_density_adjustment = a->geometry_density_adjustment;
    next.remap_error =
        fmax(marker_error,
             fmax(a->remap_error, fmax(next.remap_error, fmax(particle_error, momentum_error))));
    next.tag_elapsed = a->tag_elapsed;
    next.unresolved_seed = a->unresolved_seed;
    next.solves = a->solves;
    next.max_residual = a->max_residual;
    next.wall0 = a->wall0;
    next.rhs_calls = a->rhs_calls;
    next.tag_checks = a->tag_checks;
    next.chemistry_calls = a->chemistry_calls;
    next.timestep_attempts = a->timestep_attempts;
    next.cfl_adjustments = a->cfl_adjustments;
    next.rate_probe_calls = a->rate_probe_calls;
    next.t_exchange += a->t_exchange;
    next.t_prepare += a->t_prepare;
    next.t_coeff = a->t_coeff;
    next.t_potential = a->t_potential;
    next.t_transport = a->t_transport;
    double jump = 0;
    if (a->c.amr_case != 1) {
        a3_potential_create(&next);
        a3_potential_seed(&next, a->phi);
        a3_prepare(&next, N, M);
        a3_potential(&next);
        double scale = 1;
        for (int v = 0; v < (a->fx + 1) * (a->fy + 1); v++) {
            scale = fmax(scale, fabs(a->phi[v]));
            jump = fmax(jump, fabs(next.phi[v] - a->phi[v]));
        }
        jump /= scale;
    }
    next.max_potential_remap_change = fmax(a->max_potential_remap_change, jump);
    next.mesh_epoch = a->mesh_epoch + 1;
    next.output_mesh_epoch = a->output_mesh_epoch;
    next.regrids = a->regrids + 1;
    next.refinements = a->refinements + refine;
    next.coarsenings = a->coarsenings + coarsen;
    next.regrid_log_started = 1;
    double elapsed = wall3() - begin;
    next.t_regrid = a->t_regrid + elapsed;
    if (!a->rank) {
        char path[1200];
        snprintf(path, sizeof(path), "%s/regrid.csv", a->c.out);
        FILE *f = fopen(path, a->regrid_log_started ? "a" : "w");
        if (!f)
            a3_fail(a, "regrid audit open");
        if (!a->regrid_log_started)
            fputs(
                "time_s,step,mesh_epoch,old_leaves,new_leaves,refined_roots,coarsened_roots,ownership_changed,exact_initial,particle_error,momentum_error,potential_relative_change,wall_s\n",
                f);
        fprintf(f, "%.17g,%d,%d,%d,%d,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.9g\n", time, step,
                next.mesh_epoch, a->nleaf, next.nleaf, refine, coarsen, migration, exact_initial,
                particle_error, momentum_error, jump, elapsed);
        fclose(f);
        printf(
            "Regrid t=%.3f epoch=%d leaves=%d -> %d refine=%d coarsen=%d migrate=%d conservation=%.2g phi_change=%.2g wall=%.3fs\n",
            time, next.mesh_epoch, a->nleaf, next.nleaf, refine, coarsen, migration,
            fmax(particle_error, momentum_error), jump, elapsed);
        fflush(stdout);
    }
    free(*n);
    free(*m);
    free(a->marker);
    a3_free(a);
    *a = next;
    a3_rebase(a);
    *n = N;
    *m = M;
    return 1;
}
