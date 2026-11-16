/*
 * Experimental isothermal, night-only multi-ion closure. Quasineutral Ne equals O+ + NO+ + O2+.
 * Chemistry uses SI densities internally; stored ion densities are normalized by 1e12 m^-3.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Experimental night-only reduced multi-ion closure. NO+/O2+ have conservative
 * transport and parallel momentum; electrons use the common ion density.
 * Molecular reservoirs use the fixed canonical geometry, independently of
 * O+ AMR changes. This is not a complete multi-fluid ionosphere model.
 */
#include "amr3_night_internal.h"
#include "neutral_grid.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.14159265358979323846
/**
 * Return cumulative O+ chemical conversion loss for the global O+ conservation budget.
 */
double a3_night_loss(const A3 *a) {
    return a->night ? a->night->loss : 0;
}
/**
 * Solve perpendicular collisional gravity mobility with signed gyrofrequency and remove the
 * parallel force component.
 */
static V3 gravity_drift(V3 gravity, V3 b, double nu, double omega) {
    gravity = vsub(gravity, vmul(b, vdot(gravity, b)));
    return vmul(vsub(vmul(gravity, nu), vmul(vcross(b, gravity), omega)),
                1 / (nu * nu + omega * omega));
}
/**
 * Exchange molecular density/momentum ghost slabs; receives and sends remain valid and disjoint
 * until MPI_Waitall.
 */
static void halo(A3 *a) {
    A3Night *q = a->night;
    int width = 2 * a->fx * a->ns, count = 0;
    MPI_Request requests[8 * a->size];
    MPI_Status statuses[8 * a->size];
    /* Receives occupy disjoint peer-owned ghost slabs; sends read our owned
     * slab only. No buffer is modified before all requests have completed. */
    for (int field = 0; field < 4; field++)
        for (int peer = 0; peer < a->size; peer++)
            if (peer != a->rank) {
                double *v = field < 2 ? q->ion[field] : q->momentum[field - 2];
                int lo =
                    a->row_first > a->partition[peer] - 2 ? a->row_first : a->partition[peer] - 2;
                int hi = a->row_last < a->partition[peer + 1] + 2 ? a->row_last
                                                                  : a->partition[peer + 1] + 2;
                int send = hi > lo ? (hi - lo) * width : 0, so = send ? lo * width : 0;
                lo = a->partition[peer] > a->row_first - 2 ? a->partition[peer] : a->row_first - 2;
                hi = a->partition[peer + 1] < a->row_last + 2 ? a->partition[peer + 1]
                                                              : a->row_last + 2;
                int receive = hi > lo ? (hi - lo) * width : 0, ro = receive ? lo * width : 0;
                if (receive)
                    MPI_Irecv(v + ro, receive, MPI_DOUBLE, peer, 180 + field, MPI_COMM_WORLD,
                              requests + count++);
                if (send)
                    MPI_Isend(v + so, send, MPI_DOUBLE, peer, 180 + field, MPI_COMM_WORLD,
                              requests + count++);
            }
    MPI_Waitall(count, requests, statuses);
}
/**
 * Collect the full canonical molecular state for saved outputs/checkpoints only; normal
 * RHS/regrid paths use halos and owned migration.
 */
void a3_night_gather(A3 *a) {
    if (!a->night)
        return;
    int *c = a3_alloc(a, 2 * a->size, sizeof(int)), *o = c + a->size;
    for (int r = 0; r < a->size; r++) {
        o[r] = 2 * a->partition[r] * a->fx * a->ns;
        c[r] = 2 * (a->partition[r + 1] - a->partition[r]) * a->fx * a->ns;
    }
    for (int s = 0; s < 2; s++) {
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, a->night->ion[s], c, o, MPI_DOUBLE,
                       MPI_COMM_WORLD);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, a->night->momentum[s], c, o, MPI_DOUBLE,
                       MPI_COMM_WORLD);
    }
    free(c);
}
/**
 * Refresh current molecular Pedersen/Hall and gravity-current responses in owned plus ghost
 * cells after halo exchange.
 */
void a3_night_refresh(A3 *a) {
    if (!a->night)
        return;
    A3Night *q = a->night;
    double begin = wall3();
    halo(a);
    int j0 = a->row_first > 2 ? 2 * (a->row_first - 2) : 0;
    int j1 = a->row_last + 2 < a->ny ? 2 * (a->row_last + 2) : a->fy;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int j = j0; j < j1; j++)
        for (int i = 0; i < a->fx; i++)
            for (int k = 0; k < a->ns; k++) {
                int id = (j * a->fx + i) * a->ns + k;
                TubeCell *g = a->g->cell + id;
                double *r = a->response + 12 * id;
                double p = 0, h = 0, cx = 0, cy = 0;
                for (int s = 0; s < 2; s++) {
                    double n = q->ion[s][id], *v = q->coefficient[s] + 4 * id;
                    p += n * v[0];
                    h += n * v[1];
                    cx += n * v[2];
                    cy += n * v[3];
                }
                double xx = vdot(g->gx, g->gx), xy = vdot(g->gx, g->gy), yy = vdot(g->gy, g->gy);
                double skew = a->c.hall ? vdot(g->gx, vcross(g->b, g->gy)) : 0;
                r[4] = g->volume * p * xx;
                r[5] = g->volume * (p * xy + h * skew);
                r[6] = g->volume * (p * xy - h * skew);
                r[7] = g->volume * p * yy;
                r[10] = g->volume * cx;
                r[11] = g->volume * cy;
            }
    q->wall += wall3() - begin;
}
typedef struct {
    A3 *a;
    FILE *f;
    uint64_t hash;
    int write;
} IO;
/**
 * Read or write exactly the requested bytes while updating the checkpoint checksum; short
 * transfers are fatal.
 */
static void bytes(IO *s, void *p, size_t count) {
    if ((s->write ? fwrite(p, 1, count, s->f) : fread(p, 1, count, s->f)) != count)
        a3_fail(s->a, "night state transfer size");
    for (size_t i = 0; i < count; i++) {
        s->hash ^= ((unsigned char *)p)[i];
        s->hash *= UINT64_C(1099511628211);
    }
}
/**
 * Encode/decode an unsigned 64-bit integer in little-endian order and update the stream
 * checksum.
 */
static uint64_t integer(IO *s, uint64_t x) {
    unsigned char b[8];
    for (int k = 0; k < 8; k++)
        b[k] = (unsigned char)(x >> (8 * k));
    bytes(s, b, 8);
    x = 0;
    for (int k = 0; k < 8; k++)
        x |= (uint64_t)b[k] << (8 * k);
    return x;
}
/**
 * Transfer floating data in portable little-endian order, retaining the checksum over the actual
 * encoded bytes.
 */
static void reals(IO *s, double *p, size_t count) {
    unsigned short one = 1;
    if (*(unsigned char *)&one)
        bytes(s, p, count * 8);
    else
        for (size_t k = 0; k < count; k++) {
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
static void finish(IO *s) {
    uint64_t h = s->hash;
    if (integer(s, h) != h)
        a3_fail(s->a, "night state checksum");
    if (!s->write && (fgetc(s->f) != EOF || ferror(s->f)))
        a3_fail(s->a, "night state trailing data");
    if (fclose(s->f))
        a3_fail(s->a, "night state close");
}
/**
 * Validate and read the molecular checkpoint, numerical-source identity and binding to the
 * native checkpoint checksum.
 */
static void read_restart(A3 *a, double time) {
    A3Night *q = a->night;
    double info[9] = {0};
    if (!a->rank) {
        Solver3 check = {.rank = 0};
        IO s = {a, fopen(a->c.night_restart, "rb"), UINT64_C(14695981039346656037), 0};
        if (!s.f)
            a3_fail(a, "night restart open");
        char magic[8];
        bytes(&s, magic, 8);
        uint64_t h[7];
        for (int i = 0; i < 7; i++)
            h[i] = integer(&s, 0);
        char code[64];
        bytes(&s, code, 64);
        reals(&s, info, 9);
        if (memcmp(code, ESF_NIGHT_ID, 64))
            a3_fail(a, "night restart numerical source identity mismatch");
        if (memcmp(magic, "ESFNIT02", 8) || h[0] != (uint64_t)a->fx || h[1] != (uint64_t)a->fy ||
            h[2] != (uint64_t)a->ns || h[3] != q->count || h[4] != q->profile_hash ||
            h[5] != solver3_file_hash(&check, a->c.restart) || h[6] != 2 || info[0] != time ||
            info[1] != a->c.utc_start || info[2] != a->c.chem_Ti || info[3] != a->c.chem_Te)
            a3_fail(a, "night restart grid/profile/paired checkpoint/clock mismatch");
        for (int s0 = 0; s0 < 2; s0++) {
            reals(&s, q->ion[s0], q->count);
            for (size_t i = 0; i < q->count; i++)
                if (!isfinite(q->ion[s0][i]) || q->ion[s0][i] < 0)
                    a3_fail(a, "invalid molecular restart");
        }
        for (int s0 = 0; s0 < 2; s0++) {
            reals(&s, q->momentum[s0], q->count);
            for (size_t i = 0; i < q->count; i++)
                if (!isfinite(q->momentum[s0][i]) ||
                    (q->ion[s0][i] == 0 && q->momentum[s0][i] != 0))
                    a3_fail(a, "invalid molecular momentum restart");
        }
        finish(&s);
    }
    MPI_Bcast(info, 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int s = 0; s < 2; s++) {
        MPI_Bcast(q->ion[s], (int)q->count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(q->momentum[s], (int)q->count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    for (int i = 0; i < 8; i++)
        if (!isfinite(info[i]) || info[i] < 0)
            a3_fail(a, "invalid night restart history");
    q->loss = info[4];
    q->recomb = info[5];
    q->mass0 = info[6];
    q->max_charge_error = info[7];
    if (!isfinite(info[8]))
        a3_fail(a, "invalid molecular boundary history");
    q->boundary_total = info[8];
}
/**
 * Load frozen neutral/rate inputs, initialize NO+/O2+ reservoirs and transport coefficients, or
 * restore a paired molecular checkpoint; enforce the night-only applicability guard.
 */
void a3_night_create(A3 *a, double time) {
    if (!a->c.night_chemistry)
        return;
    A3Night *q = a3_alloc(a, 1, sizeof(*q));
    a->night = q;
    q->count = (size_t)a->fx * a->fy * a->ns;
    Solver3 check = {.rank = a->rank};
    q->profile_hash = solver3_file_hash(&check, a->c.neutral_profile);
    q->rate = a3_alloc(a, 2 * q->count, sizeof(double));
    q->electron = a3_alloc(a, q->count, sizeof(double));
    for (int s = 0; s < 2; s++) {
        q->ion[s] = a3_alloc(a, q->count, sizeof(double));
        q->momentum[s] = a3_alloc(a, q->count, sizeof(double));
        q->nu[s] = a3_alloc(a, q->count, sizeof(double));
        q->coefficient[s] = a3_alloc(a, 4 * q->count, sizeof(double));
    }
    double *spatial = NULL;
    int spatial_status =
        neutral_grid_read(a->c.neutral_profile, a->g, a->c.chem_Ti, a->c.chem_Te, &spatial);
    if (spatial_status < 0 || !!spatial_status != !!a->c.spatial_neutral)
        a3_fail(a, "neutral grid input or mode changed");
    double table[4096][5];
    int nt = 0;
    char line[2048];
    FILE *f = NULL;
    if (!spatial) {
        f = fopen(a->c.neutral_profile, "r");
        if (!f)
            a3_fail(a, "neutral table open");
    }
    while (f && fgets(line, sizeof(line), f)) {
        if (!strncmp(line, "height_m,", 9) || line[0] == '#' || line[0] == '\n')
            continue;
        if (nt >= 4096 || sscanf(line, "%lf,%lf,%lf,%lf,%lf", table[nt], table[nt] + 1,
                                 table[nt] + 2, table[nt] + 3, table[nt] + 4) != 5)
            a3_fail(a, "neutral table columns");
        for (int j = 0; j < 5; j++)
            if (!isfinite(table[nt][j]) || table[nt][j] <= 0)
                a3_fail(a, "neutral table positive values");
        if (nt && table[nt][0] <= table[nt - 1][0])
            a3_fail(a, "neutral table height order");
        nt++;
    }
    if (f)
        fclose(f);
    if (!spatial && nt < 2)
        a3_fail(a, "neutral table empty");
    const double qe = 1.602176634e-19, amu = 1.66053906660e-27, me = 9.1093837139e-31;
    for (size_t id = 0; id < q->count; id++) {
        TubeCell *g = a->g->cell + id;
        double h = g->height, v[4];
        if (spatial) {
            for (int j = 0; j < 4; j++)
                v[j] = spatial[6 * id + j];
        } else {
            if (h < table[0][0] || h > table[nt - 1][0])
                a3_fail(a, "neutral table height coverage");
            int lo = 0, hi = nt - 1;
            while (hi - lo > 1) {
                int mid = (lo + hi) / 2;
                if (table[mid][0] > h)
                    hi = mid;
                else
                    lo = mid;
            }
            double w = (h - table[lo][0]) / (table[hi][0] - table[lo][0]);
            for (int j = 0; j < 4; j++)
                v[j] = exp((1 - w) * log(table[lo][j + 1]) + w * log(table[hi][j + 1]));
        }
        double O = v[0] * 1e-6, N2 = v[1] * 1e-6, O2 = v[2] * 1e-6, Tn = v[3], Ti = a->c.chem_Ti,
               Te = a->c.chem_Te;
        NightRates r;
        if (night_rates(Ti, Te, Tn, v[1], v[2], &r))
            a3_fail(a, "reaction fit temperature range");
        q->rate[2 * id] = r.to_no;
        q->rate[2 * id + 1] = r.to_o2;
        q->alpha[0] = r.alpha_no;
        q->alpha[1] = r.alpha_o2;
        double Tr = .5 * (Ti + Tn),
               nu[2] = {4.27e-10 * O2 + 2.44e-10 * O + 4.34e-10 * N2,
                        2.59e-11 * O2 * sqrt(Tr) * pow(1 - .073 * log10(Tr), 2) + 2.31e-10 * O +
                            4.13e-10 * N2};
        double nue = 4 * (2.33e-11 * N2 * Te * (1 - 1.21e-4 * Te) +
                          1.82e-10 * O2 * sqrt(Te) * (1 + .036 * sqrt(Te)) +
                          8.9e-11 * O * sqrt(Te) * (1 + 5.7e-4 * Te));
        double oe = qe * g->B / me, ke = oe / nue, pe = ke / (1 + ke * ke);
        V3 grav = vmul(vunit(g->position), -9.8 * pow(6371200. / vnorm(g->position), 2));
        for (int s = 0; s < 2; s++) {
            q->nu[s][id] = nu[s];
            double om = qe * g->B / ((s ? 32. : 30.) * amu), ki = om / nu[s],
                   *c = q->coefficient[s] + 4 * id;
            c[0] = qe / g->B * 1e12 * (ki / (1 + ki * ki) + pe);
            c[1] = qe / g->B * 1e12 * (1 / (1 + ki * ki) - 1 / (1 + ke * ke));
            V3 j = vmul(
                vsub(gravity_drift(grav, g->b, nu[s], om), gravity_drift(grav, g->b, nue, -oe)),
                qe * 1e12);
            c[2] = vdot(j, g->gx);
            c[3] = vdot(j, g->gy);
        }
        q->ion[0][id] = g->nE / 1e12;
        q->mass0 += g->volume * q->ion[0][id];
    }
    free(spatial);
    /* Require all modeled 90–150km molecular-E cells to remain in geometric
     * Earth shadow. No arbitrary sunlight switch or fake photo-production. */
    for (double t = time;; t = fmin(t + 1800, a->c.end)) {
        V3 sun = {night_solar_cosine(a->c.utc_start + t, 0, 0),
                  night_solar_cosine(a->c.utc_start + t, 0, PI / 2),
                  night_solar_cosine(a->c.utc_start + t, PI / 2, 0)};
        for (size_t i = 0; i < q->count; i++) {
            TubeCell *g = a->g->cell + i;
            if (g->height > 150e3)
                continue;
            double radius = vnorm(g->position), mu = vdot(g->position, sun) / radius;
            if (mu >= 0 || radius * radius * (1 - mu * mu) > 6371200. * 6371200.)
                a3_fail(
                    a, "night-only E-region closure requested during modeled E-layer illumination");
        }
        if (t >= a->c.end)
            break;
    }
    if (a->c.restart[0])
        read_restart(a, time);
    a3_night_refresh(a);
}
/**
 * Back up molecular state and reaction budgets before an attempted coupled step.
 */
void a3_night_begin(A3 *a) {
    if (!a->night)
        return;
    A3Night *q = a->night;
    size_t count = night_owned(a), off = night_start(a);
    a3_molecular_workspace(a);
    for (int s = 0; s < 2; s++) {
        memcpy(q->backup[s], q->ion[s] + off, count * sizeof(double));
        memcpy(q->backup_m[s], q->momentum[s] + off, count * sizeof(double));
    }
    q->old_loss = q->loss;
    q->old_recomb = q->recomb;
}
/**
 * Restore the pre-attempt molecular state and reaction budgets when a CFL or positivity check
 * requires another attempt.
 */
void a3_night_restore(A3 *a) {
    if (!a->night)
        return;
    A3Night *q = a->night;
    for (int s = 0; s < 2; s++) {
        memcpy(q->ion[s] + night_start(a), q->backup[s], night_owned(a) * sizeof(double));
        memcpy(q->momentum[s] + night_start(a), q->backup_m[s], night_owned(a) * sizeof(double));
    }
    q->loss = q->old_loss;
    q->recomb = q->old_recomb;
}
/**
 * Apply local reactions to disjoint canonical children of actual O+ volumes, then volume-average
 * O+ back to its owner. Thread-local budgets are reduced; fatal handling stays outside workers.
 */
static void macro_chemistry(A3 *a, double *n, double *m, double dt, double budget[2],
                            double *error) {
    A3Night *q = a->night;
    double converted = 0, recombined = 0, worst = 0;
    int invalid = 0;
/* Each active finite volume maps to disjoint canonical children. MPI
 * reductions and fatal error handling stay outside worker threads. */
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) reduction(+ : converted, recombined)         \
    reduction(max : worst) reduction(| : invalid)
#endif
    for (int slot = 0; slot < a->owned_count; slot++)
        for (int k = 0; k < a->ns; k++) {
            int p = a->global[slot];
            if (!a3_active(a, p, k))
                continue;
            A3Leaf l = a3_density_leaf(a, p, k);
            int scale = a3_scale(l.level), t = slot * a->ns + k;
            double initial = n[t], mass = 0;
            for (int v = 0; v < scale; v++)
                for (int u = 0; u < scale; u++) {
                    int id = tube_cell_id(a->g, scale * l.i + u, scale * l.j + v, k);
                    double o = initial * 1e12, no = q->ion[0][id] * 1e12, o2 = q->ion[1][id] * 1e12,
                           total = o + no + o2;
                    NightRates rates = {q->rate[2 * id], q->rate[2 * id + 1], q->alpha[0],
                                        q->alpha[1]};
                    NightBudget b;
                    if (night_react(&o, &no, &o2, rates, dt, &b)) {
                        invalid = 1;
                        continue;
                    }
                    worst =
                        fmax(worst, fabs((o + no + o2 + b.recombined - total) / fmax(total, 1.)));
                    double V = a->g->cell[id].volume;
                    mass += V * o / 1e12;
                    converted += b.converted / 1e12 * V;
                    recombined += b.recombined / 1e12 * V;
                    q->momentum[0][id] *= q->ion[0][id] > 0 ? (no / 1e12) / q->ion[0][id] : 0;
                    q->momentum[1][id] *= q->ion[1][id] > 0 ? (o2 / 1e12) / q->ion[1][id] : 0;
                    q->ion[0][id] = no / 1e12;
                    q->ion[1][id] = o2 / 1e12;
                }
            double next = mass / a->metric[p * a->ns + k].volume,
                   factor = initial > 0 ? next / initial : 0;
            n[t] = next;
            m[t] *= factor;
            if (a->marker)
                a->marker[t] *= factor;
        }
    if (invalid)
        a3_fail(a, "invalid multilevel chemistry state");
    budget[0] += converted;
    budget[1] += recombined;
    *error = fmax(*error, worst);
}
/**
 * Execute one chemistry split operator, check its charge budget, reduce global
 * conversion/recombination totals. The driver defers molecular halo/response
 * refresh to a3_prepare, immediately before the next current-state consumer.
 */
void a3_night_step(A3 *a, double *n, double *m, double dt) {
    if (!a->night || !a->c.night_reactions)
        return;
    A3Night *q = a->night;
    double begin = wall3(), budget[2] = {0}, error = 0;
    a->chemistry_calls++;
    if (a->c.amr_levels == 3)
        macro_chemistry(a, n, m, dt, budget, &error);
    else {
        for (int j = a->row_first; j < a->row_last; j++)
            for (int i = 0; i < a->nx; i++)
                for (int k = 0; k < a->ns; k++) {
                    int root = j * a->nx + i, refined = a->density_refined[root * a->ns + k],
                        primary = a->first[root];
                    double mass = 0, volume = 0, old = 0;
                    for (int d = 0; d < 4; d++) {
                        int canonical = ((2 * j + (d >> 1)) * a->fx + 2 * i + (d & 1)) * a->ns + k;
                        int leaf = refined ? primary + d : primary,
                            slot = a->slot[leaf] * a->ns + k;
                        double initial = n[slot], o = initial * 1e12,
                               no = q->ion[0][canonical] * 1e12, o2 = q->ion[1][canonical] * 1e12,
                               total = o + no + o2;
                        NightRates rates = {q->rate[2 * canonical], q->rate[2 * canonical + 1],
                                            q->alpha[0], q->alpha[1]};
                        NightBudget b;
                        if (night_react(&o, &no, &o2, rates, dt, &b))
                            a3_fail(a, "invalid chemistry state");
                        error = fmax(error,
                                     fabs((o + no + o2 + b.recombined - total) / fmax(total, 1.)));
                        double V = a->g->cell[canonical].volume, next = o / 1e12;
                        /* Density-only local chemistry holds each species velocity. The
                         * fixed neutral bath supplies/removes reaction momentum; no claim
                         * of reactive momentum conservation between charged species. */
                        q->momentum[0][canonical] *=
                            q->ion[0][canonical] > 0 ? (no / 1e12) / q->ion[0][canonical] : 0;
                        q->momentum[1][canonical] *=
                            q->ion[1][canonical] > 0 ? (o2 / 1e12) / q->ion[1][canonical] : 0;
                        q->ion[0][canonical] = no / 1e12;
                        q->ion[1][canonical] = o2 / 1e12;
                        budget[0] += b.converted / 1e12 * V;
                        budget[1] += b.recombined / 1e12 * V;
                        if (refined) {
                            double factor = initial > 0 ? next / initial : 0;
                            n[slot] = next;
                            m[slot] *= factor;
                            if (a->marker)
                                a->marker[slot] *= factor;
                        } else {
                            mass += next * V;
                            volume += V;
                            old = initial;
                        }
                    }
                    if (!refined) {
                        int slot = a->slot[primary] * a->ns + k;
                        double next = mass / volume, factor = old > 0 ? next / old : 0;
                        n[slot] = next;
                        m[slot] *= factor;
                        if (a->marker)
                            a->marker[slot] *= factor;
                    }
                }
    }
    MPI_Allreduce(MPI_IN_PLACE, budget, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (error > 5e-14)
        a3_fail(a, "reaction charge budget failed");
    q->max_charge_error = fmax(q->max_charge_error, error);
    q->loss += budget[0];
    q->recomb += budget[1];
    q->wall += wall3() - begin;
    /* Only owned species changed. Both callers return to the driver, whose
     * next current-state consumer runs a3_prepare before potential/transport
     * or tagging; output also starts with a full RHS. Prepare exchanges the
     * new molecular halo and refreshes every conductivity response. Doing
     * that exchange/refresh here as well repeats identical work before any
     * response is read. No old response may be consumed between this return
     * and prepare. Chemistry/rollback need only owned molecular state. */
}
/**
 * Write a molecular checkpoint bound to the just-written native checkpoint, molecular frames,
 * metadata and total-ion reaction/boundary budgets.
 */
void a3_night_snapshot(A3 *a, const double *n, double time, int frame, double mass0) {
    if (!a->night)
        return;
    A3Night *q = a->night;
    a3_night_gather(a);
    double omass = a3_mass(a, n);
    if (a->rank)
        return;
    double mmass = 0;
    for (size_t k = 0; k < q->count; k++)
        mmass += (q->ion[0][k] + q->ion[1][k]) * a->g->cell[k].volume;
    char path[1200], part[1200], native[1200];
    snprintf(path, sizeof(path), "%s/night_state_latest.bin", a->c.out);
    snprintf(part, sizeof(part), "%s/night_state_latest.part", a->c.out);
    snprintf(native, sizeof(native), "%s/checkpoint_latest.bin", a->c.out);
    Solver3 check = {.rank = 0};
    IO s = {a, fopen(part, "wb"), UINT64_C(14695981039346656037), 1};
    if (!s.f)
        a3_fail(a, "night checkpoint open");
    char magic[] = "ESFNIT02";
    bytes(&s, magic, 8);
    uint64_t h[] = {
        a->fx, a->fy, a->ns, q->count, q->profile_hash, solver3_file_hash(&check, native), 2};
    for (int k = 0; k < 7; k++)
        integer(&s, h[k]);
    char code[] = ESF_NIGHT_ID;
    bytes(&s, code, 64);
    double info[] = {time,      a->c.utc_start, a->c.chem_Ti,        a->c.chem_Te,     q->loss,
                     q->recomb, q->mass0,       q->max_charge_error, q->boundary_total};
    reals(&s, info, 9);
    for (int k = 0; k < 2; k++)
        reals(&s, q->ion[k], q->count);
    for (int k = 0; k < 2; k++)
        reals(&s, q->momentum[k], q->count);
    finish(&s);
    if (rename(part, path))
        a3_fail(a, "night checkpoint rename");
    snprintf(path, sizeof(path), "%s/molecular_%04d.bin", a->c.out, frame);
    FILE *f = fopen(path, "wb");
    if (!f)
        a3_fail(a, "molecular output open");
    for (int k = 0; k < 2; k++)
        if (fwrite(q->ion[k], sizeof(double), q->count, f) != q->count)
            a3_fail(a, "molecular output write");
    fclose(f);
    snprintf(path, sizeof(path), "%s/molecular_momentum_%04d.bin", a->c.out, frame);
    f = fopen(path, "wb");
    if (!f)
        a3_fail(a, "molecular momentum output open");
    for (int k = 0; k < 2; k++)
        if (fwrite(q->momentum[k], sizeof(double), q->count, f) != q->count)
            a3_fail(a, "molecular momentum output write");
    fclose(f);
    snprintf(path, sizeof(path), "%s/night_diagnostics.csv", a->c.out);
    f = fopen(path, frame ? "a" : "w");
    if (!f)
        a3_fail(a, "night diagnostics open");
    if (!frame)
        fputs(
            "time_s,utc_unix_s,Oplus_particle_loss,recombined_electrons,molecular_particle_mass,total_ion_relative_balance,maximum_kernel_charge_error,chemistry_wall_s,molecular_boundary_particles,molecular_transport_wall_s\n",
            f);
    double balance =
        (omass + mmass - mass0 - q->mass0 - a->boundary_total - q->boundary_total + q->recomb) /
        (mass0 + q->mass0);
    fprintf(f, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.9g,%.17g,%.9g\n", time,
            a->c.utc_start + time, q->loss, q->recomb, mmass, balance, q->max_charge_error, q->wall,
            q->boundary_total, q->transport_wall);
    fclose(f);
    if (fabs(balance) > 1e-10)
        a3_fail(a, "global ion/electron reaction budget failed");
    snprintf(path, sizeof(path), "%s/night_metadata.json", a->c.out);
    f = fopen(path, "w");
    if (!f)
        a3_fail(a, "night metadata open");
    fprintf(
        f,
        "{\"night_model\":\"experimental_transported_multi_ion_v2\",\"utc_start_unix_s\":%.17g,\"Ti_K\":%.17g,\"Te_K\":%.17g,\"molecular_format\":\"float64 [2,fine_ny,fine_nx,ns], NO+ then O2+, normalized1e12m^-3\",\"initial_molecular_composition\":\"all prescribed E population initially NO+\",\"molecular_transport\":true,\"common_electron_pressure\":true,\"molecular_perpendicular_scheme\":\"donor_cell_on_canonical_fine_grid\",\"perpendicular_molecular_pressure_diffusion\":false,\"reactive_momentum\":\"species_velocity_held_in_fixed_neutral_bath\",\"night_algorithm\":2,\"photoionization\":false,\"night_guard\":\"modeled90-150km E layer in geometric Earth shadow, checked every1800s\",\"clock_is_MLT\":false,\"Oplus_chemical_loss\":%.17g,\"recombined_electrons\":%.17g,\"total_ion_relative_balance\":%.17g,\"maximum_kernel_charge_error\":%.17g,\"chemistry_wall_s\":%.9g}\n",
        a->c.utc_start, a->c.chem_Ti, a->c.chem_Te, q->loss, q->recomb, balance,
        q->max_charge_error, q->wall);
    if (a->c.spatial_neutral) {
        /* Separate metadata keeps the legacy schema readable. */
        FILE *bf;
        char bp[1200];
        snprintf(bp, sizeof(bp), "%s/neutral_grid_metadata.json", a->c.out);
        bf = fopen(bp, "w");
        if (!bf)
            a3_fail(a, "spatial neutral metadata");
        fprintf(
            bf,
            "{\"spatial_neutral\":true,\"profile_hash_fnv64\":\"%016llx\",\"frozen_at_start_UTC\":true,\"initial_plasma_still_idealized\":true,\"Ti_K\":%.17g,\"Te_K\":%.17g}\n",
            (unsigned long long)a->c.neutral_grid_hash, a->c.chem_Ti, a->c.chem_Te);
        fclose(bf);
    }
    fclose(f);
}
/**
 * Release all molecular full-grid and owned-work arrays, including stage-local ETD/log caches.
 */
void a3_night_free(A3 *a) {
    if (a->night)
        for (int f = 0; f < 3; f++) {
            free(a->night->electron_average[f]);
            free(a->night->electron_ready[f]);
        }
    if (a->night)
        for (int s = 0; s < 2; s++)
            for (int f = 0; f < 3; f++)
                free(a->night->etd_weight[s][f]);
    if (!a->night)
        return;
    A3Night *q = a->night;
    for (int s = 0; s < 2; s++) {
        free(q->ion[s]);
        free(q->momentum[s]);
        free(q->nu[s]);
        free(q->coefficient[s]);
        free(q->backup[s]);
        free(q->backup_m[s]);
        free(q->old_n[s]);
        free(q->old_m[s]);
        free(q->dn[s]);
        free(q->dm[s]);
        free(q->first_m[s]);
        free(q->rates[s]);
    }
    for (int t = 0; t < 3; t++) {
        free(q->face[t]);
        free(q->transport_log[t]);
    }
    free(q->electron);
    free(q->rate);
    free(q);
    a->night = NULL;
}
