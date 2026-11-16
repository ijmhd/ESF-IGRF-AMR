/*
 * NO+ and O2+ transport on canonical finest cells, with separate parallel momenta and common
 * electron pressure. Canonical geometry is fixed even when O+ adapts.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Canonical-fine finite volumes for NO+/O2+, separate parallel momenta.
 * Species mobility across B; common Ne electron pressure along B.
 * Perpendicular molecular pressure diffusion and reactive momentum exchange
 * are outside this reduced, fixed-neutral-bath/isothermal model.
 */
#include "amr3_night_internal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
static const double amu = 1.66053906660e-27, kb = 1.380649e-23, qe = 1.602176634e-19;
/**
 * Clamp a cell index to the nearest valid boundary cell.
 */
static int clamp(int x, int n) {
    return x < 0 ? 0 : x >= n ? n - 1 : x;
}
/**
 * Flatten a canonical molecular cell with periodic-x or clamped physical indices.
 */
static int id(A3 *a, int i, int j, int k) {
    if (a->c.periodic_x) {
        i %= a->fx;
        if (i < 0)
            i += a->fx;
    }
    return (clamp(j, a->fy) * a->fx + clamp(i, a->fx)) * a->ns + clamp(k, a->ns);
}
/**
 * Flatten node coordinates with the along-field coordinate contiguous.
 */
static int nid(A3 *a, int i, int j, int k) {
    return (j * (a->fx + 1) + i) * (a->ns + 1) + k;
}
/**
 * Average the four geometric vertices of a canonical face in ECEF metres.
 */
static V3 center(A3 *a, int axis, int p[3]) {
    int b = (axis + 1) % 3, c = (axis + 2) % 3;
    V3 x = {0, 0, 0};
    for (int u = 0; u < 2; u++)
        for (int v = 0; v < 2; v++) {
            int q[3] = {p[0], p[1], p[2]};
            q[b] += u;
            q[c] += v;
            x = vadd(x, a->g->node[nid(a, q[0], q[1], q[2])]);
        }
    return vmul(x, .25);
}
/**
 * Build owned canonical face mobility/gravity coefficients, and rebuild only when the owned slab
 * changes.
 */
void a3_molecular_geometry(A3 *a) {
    A3Night *q = a->night;
    if (!q || (q->face[0] && q->face_first == a->row_first && q->face_last == a->row_last))
        return;
    Igrf field;
    if (solver3_load_field(&a->c, &field))
        a3_fail(a, "molecular magnetic field");
    int n[3] = {a->fx, a->fy, a->ns};
    for (int axis = 0; axis < 3; axis++) {
        free(q->face[axis]);
        int ni = a->fx + (axis == 0), nk = a->ns + (axis == 2),
            nj = 2 * (a->row_last - a->row_first) + (axis == 1);
        q->face[axis] = a3_alloc(a, (size_t)ni * nk * nj * 6, sizeof(double));
        for (int j = 2 * a->row_first; j < 2 * a->row_last + (axis == 1); j++)
            for (int i = 0; i < ni; i++)
                for (int k = 0; k < nk; k++) {
                    int p[3] = {i, j, k}, b = (axis + 1) % 3, c = (axis + 2) % 3;
                    V3 e[3] = {{0, 0, 0}}, x = center(a, axis, p);
                    for (int t = 0; t < 2; t++) {
                        int l[3] = {i, j, k}, r[3] = {i, j, k};
                        l[c] += t;
                        r[c] += t;
                        r[b]++;
                        e[b] = vadd(e[b], vmul(vsub(a->g->node[nid(a, r[0], r[1], r[2])],
                                                    a->g->node[nid(a, l[0], l[1], l[2])]),
                                               .5 * n[b]));
                        l[0] = r[0] = i;
                        l[1] = r[1] = j;
                        l[2] = r[2] = k;
                        l[b] += t;
                        r[b] += t;
                        r[c]++;
                        e[c] = vadd(e[c], vmul(vsub(a->g->node[nid(a, r[0], r[1], r[2])],
                                                    a->g->node[nid(a, l[0], l[1], l[2])]),
                                               .5 * n[c]));
                    }
                    int pp[3] = {i, j, k}, pm[3] = {i, j, k};
                    pp[axis] = p[axis] < n[axis] ? p[axis] + 1 : p[axis];
                    pm[axis] = p[axis] > 0 ? p[axis] - 1 : p[axis];
                    e[axis] = vmul(vsub(center(a, axis, pp), center(a, axis, pm)),
                                   (double)n[axis] / (pp[axis] - pm[axis]));
                    double det = vdot(e[0], vcross(e[1], e[2]));
                    if (!isfinite(det) || det == 0)
                        a3_fail(a, "molecular face Jacobian");
                    V3 gx = vmul(vcross(e[1], e[2]), 1 / det),
                       gy = vmul(vcross(e[2], e[0]), 1 / det);
                    V3 S = vmul(vcross(e[b], e[c]), (det > 0 ? 1. : -1.) / (n[b] * n[c]));
                    V3 bf = igrf_field(&field, x), bh = vunit(e[2]);
                    double B = vnorm(bf);
                    if (vdot(bh, bf) < 0)
                        bh = vmul(bh, -1);
                    gx = vsub(gx, vmul(bh, vdot(gx, bh)));
                    gy = vsub(gy, vmul(bh, vdot(gy, bh)));
                    V3 grav = vmul(vunit(x), -9.8 * pow(6371200. / vnorm(x), 2));
                    grav = vsub(grav, vmul(bh, vdot(grav, bh)));
                    int L = id(a, i - (axis == 0), j - (axis == 1), k - (axis == 2)),
                        R = id(a, i, j, k);
                    TubeFace *f = a->g->face[axis] + tube_face_id(a->g, axis, i, j, k);
                    double *out = q->face[axis] + 6 * (((j - 2 * a->row_first) * ni + i) * nk + k);
                    for (int s = 0; s < 2; s++) {
                        double nu = .5 * (q->nu[s][L] + q->nu[s][R]),
                               om = qe * B / ((s ? 32 : 30) * amu), h, ped;
                        night_mobility(nu, om, &h, &ped);
                        V3 drift = vmul(vsub(vmul(grav, nu), vmul(vcross(bh, grav), om)),
                                        1 / (nu * nu + om * om));
                        out[3 * s] = h * f->e[0] - ped / B * vdot(S, gx);
                        out[3 * s + 1] = h * f->e[1] - ped / B * vdot(S, gy);
                        out[3 * s + 2] = vdot(S, drift);
                    }
                }
    }
    q->face_first = a->row_first;
    q->face_last = a->row_last;
}
/**
 * Resize owned molecular stage, backup, rate, pressure-log, averaged-Ne and same-attempt ETD
 * work arrays after an ownership change. No global-sized cache is introduced.
 */
void a3_molecular_workspace(A3 *a) {
    A3Night *q = a->night;
    size_t count = night_owned(a);
    if (q->back_count == count)
        return;
    for (int s = 0; s < 2; s++) {
        double **list[] = {&q->backup[s], &q->backup_m[s], &q->old_n[s],   &q->old_m[s],
                           &q->dn[s],     &q->dm[s],       &q->first_m[s], &q->rates[s]};
        for (int t = 0; t < 8; t++) {
            free(*list[t]);
            *list[t] = a3_alloc(a, count, sizeof(double));
        }
    }
    for (int f = 0; f < 3; f++) {
        free(q->transport_log[f]);
        q->transport_log[f] = a3_alloc(a, count, sizeof(double));
    }
    for (int s = 0; s < 2; s++)
        for (int f = 0; f < 3; f++) {
            free(q->etd_weight[s][f]);
            q->etd_weight[s][f] = a3_alloc(a, count, sizeof(double));
        }
    for (int level = a->c.amr_levels == 3 ? -1 : 0; level <= 1; level++) {
        int index = level + 1, scale = a3_scale(level);
        size_t cells = count / (scale * scale);
        free(q->electron_average[index]);
        free(q->electron_ready[index]);
        q->electron_average[index] = a3_alloc(a, cells, sizeof(double));
        q->electron_ready[index] = a3_alloc(a, cells, 1);
        q->electron_count[index] = cells;
    }
    q->back_count = count;
}
/**
 * Reconstruct canonical O+ and add both molecular ions to form current-stage quasineutral
 * electron density and invalidate every averaged-Ne sample. Sampling has lazy caches; this
 * loop stays on the main thread.
 */
void a3_night_electrons(A3 *a) {
    if (!a->night)
        return;
    A3Night *q = a->night;
    a3_molecular_workspace(a);
    /* This is the state-version boundary for every averaged-Ne sample.
     * Never keep these validity bits across a chemistry/predictor RHS. */
    for (int level = a->c.amr_levels == 3 ? -1 : 0; level <= 1; level++)
        memset(q->electron_ready[level + 1], 0, q->electron_count[level + 1]);
    for (int j = 2 * a->row_first; j < 2 * a->row_last; j++)
        for (int i = 0; i < a->fx; i++)
            for (int k = 0; k < a->ns; k++) {
                int t = id(a, i, j, k);
                q->electron[t] = a3_sample(a, 0, 1, i, j, k) + q->ion[0][t] + q->ion[1][t];
            }
}
/**
 * Return the volume-weighted current-stage electron density over canonical children of a
 * requested density cell. Memoize the original sum/division only for owned cells of this RHS
 * state; even the finest-level multiply/divide is retained for arithmetic parity. This lazy
 * function is called by the main-thread O+ pressure/parallel-flux path, not worker kernels.
 */
double a3_night_electron(A3 *a, int level, int i, int j, int k) {
    A3Night *q = a->night;
    int scale = a3_scale(level);
    int first = 2 * a->row_first / scale, last = 2 * a->row_last / scale, nx = a3_nx(a, level);
    int cached = level >= -1 && level <= 1 && j >= first && j < last && i >= 0 && i < nx &&
                 k >= 0 && k < a->ns;
    size_t local = cached ? ((size_t)(j - first) * nx + i) * a->ns + k : 0;
    if (cached && q->electron_ready[level + 1][local])
        return q->electron_average[level + 1][local];
    double mass = 0, volume = 0;
    for (int v = 0; v < scale; v++)
        for (int u = 0; u < scale; u++) {
            int t = id(a, scale * i + u, scale * j + v, k);
            double V = a->g->cell[t].volume;
            mass += q->electron[t] * V;
            volume += V;
        }
    double average = mass / volume;
    if (cached) {
        q->electron_average[level + 1][local] = average;
        q->electron_ready[level + 1][local] = 1;
    }
    return average;
}
/**
 * Read the cached canonical physical parallel segment length.
 */
static double length(A3 *a, int i, int j, int k) {
    return a->parallel_length[1][id(a, i, j, k)];
}
/**
 * Use current-stage cached logarithms for positive states and the original bounded secant at
 * exact vacuum.
 */
static double cached_jump(double left, double right, double log_left, double log_right) {
    /* Preserve the original exact-vacuum secant; no density floor. */
    return left > 0 && right > 0 ? log_right - log_left : night_log_jump(left, right);
}
/**
 * Compute a nonuniform centered/one-sided parallel logarithmic gradient from the current-stage
 * pressure cache.
 */
static double slope(A3 *a, const double *values, const double *logs, int i, int j, int k) {
    int lo = k ? k - 1 : k, hi = k + 1 < a->ns ? k + 1 : k, t = id(a, i, j, k), l = id(a, i, j, lo),
        r = id(a, i, j, hi);
    size_t off = night_start(a);
    double h = length(a, i, j, k), hl = .5 * (h + length(a, i, j, lo)),
           hr = .5 * (h + length(a, i, j, hi));
    double dl = cached_jump(values[l], values[t], logs[l - off], logs[t - off]);
    double dr = cached_jump(values[t], values[r], logs[t - off], logs[r - off]);
    if (k == 0)
        return dr / hr;
    if (k == a->ns - 1)
        return dl / hl;
    return (hr / hl * dl + hl / hr * dr) / (hl + hr);
}
/**
 * Return the sign-consistent slope with smallest magnitude; return zero at an extremum.
 */
static double minmod(double l, double m, double r) {
    return l * m > 0 && l * r > 0 ? copysign(fmin(fabs(l), fmin(fabs(m), fabs(r))), l) : 0;
}
/**
 * Reconstruct a bounded parallel molecular density or velocity face state using the three-point
 * slope limiter.
 */
static double state(A3 *a, int s, int i, int j, int k, int velocity, int side) {
    A3Night *q = a->night;
    double v[3];
    for (int t = -1; t <= 1; t++) {
        int x = id(a, i, j, k + t);
        double n = q->ion[s][x];
        v[t + 1] = velocity ? (n > 0 ? q->momentum[s][x] / n : 0) : n;
    }
    return v[1] + .5 * side * minmod(1.3 * (v[1] - v[0]), .5 * (v[2] - v[0]), 1.3 * (v[2] - v[1]));
}
/**
 * Evaluate both species' conservative canonical fluxes, parallel pressure/gravity sources and
 * stability rates. Threaded sweeps assign disjoint full grid lines; MPI reductions remain on the
 * main thread.
 */
void a3_night_transport(A3 *a) {
    if (!a->night)
        return;
    double begin = wall3();
    A3Night *q = a->night;
    a3_molecular_workspace(a);
    a3_molecular_geometry(a);
    size_t off = night_start(a), count = night_owned(a);
    q->boundary_rate = 0;
    double maximum_rate = 0;
    /* These are the current-stage states produced by a3_prepare. Refill after
     * every predictor/corrector state change, including after repartition. */
    if (a->c.parallel && !a->rate_only) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (size_t t = 0; t < count; t++) {
            for (int f = 0; f < 2; f++) {
                double n = q->ion[f][off + t];
                q->transport_log[f][t] = n > 0 ? log(n) : 0;
            }
            double ne = q->electron[off + t];
            q->transport_log[2][t] = ne > 0 ? log(ne) : 0;
        }
    }
    for (int s = 0; s < 2; s++) {
        if (!a->rate_only) {
            memset(q->dn[s], 0, count * sizeof(double));
            memset(q->dm[s], 0, count * sizeof(double));
        }
        memset(q->rates[s], 0, count * sizeof(double));
    }
    for (int axis = 0; axis < 3; axis++) {
        int ni = a->fx + (axis == 0), nk = a->ns + (axis == 2);
        int imax = a->fx + (axis == 0 && !a->c.periodic_x);
        /* A thread owns all x-lines in one j row, all y-lines in one i
         * column, or one (j,i) field line. Different workers still write
         * disjoint cells during an axis sweep, but k stays innermost and
         * contiguous in memory. Splitting (j,k)/(i,k) into separate work
         * items would create a needlessly strided innermost face loop. */
        int j0 = 2 * a->row_first, j1 = 2 * a->row_last,
            lines = axis == 0   ? (j1 - j0)
                    : axis == 1 ? a->fx
                                : (j1 - j0) * a->fx;
        double boundary = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : boundary)
#endif
        for (int line = 0; line < lines; line++) {
            int end = axis == 0 ? imax : axis == 1 ? j1 + 1 : nk;
            int begin = axis == 1 ? j0 : 0;
            for (int face = begin; face < end; face++) {
                int j = axis == 1 ? face : axis == 0 ? j0 + line : j0 + line / a->fx;
                int i = axis == 0 ? face : axis == 1 ? line : line % a->fx;
                double *d = a->potential_gradient[axis] + 2 * (j * (a->fx + 1) + i);
                for (int k = axis == 2 ? face : 0; k < (axis == 2 ? face + 1 : nk); k++) {
                    int wall = (axis == 0 && !a->c.periodic_x && (i == 0 || i == a->fx)) ||
                               (axis == 1 && (j == 0 || j == a->fy)) ||
                               (axis == 2 && (k == 0 || k == a->ns));
                    int open = axis == 0   ? a->c.regional_x
                               : axis == 1 ? a->c.open_apex
                                           : a->c.open_ends;
                    int L = id(a, i - (axis == 0), j - (axis == 1), k - (axis == 2)),
                        R = id(a, i, j, k);
                    int left_exists = !(axis == 0 && !a->c.periodic_x && i == 0) &&
                                      !(axis == 1 && j == 0) && !(axis == 2 && k == 0);
                    int right_exists = !(axis == 0 && !a->c.periodic_x && i == a->fx) &&
                                       !(axis == 1 && j == a->fy) && !(axis == 2 && k == a->ns);
                    int ownL = left_exists && (size_t)L >= off && (size_t)L < off + count,
                        ownR = right_exists && (size_t)R >= off && (size_t)R < off + count;
                    if (!ownL && !ownR)
                        continue;
                    const double *c =
                        q->face[axis] + 6 * (((j - 2 * a->row_first) * ni + i) * nk + k);
                    TubeFace *f = a->g->face[axis] + tube_face_id(a->g, axis, i, j, k);
                    for (int s = 0; s < 2; s++) {
                        double v = c[3 * s] * d[0] + c[3 * s + 1] * d[1] + c[3 * s + 2];
                        if (wall && !open)
                            v = 0;
                        double F = 0, H = 0, rate = fabs(v);
                        if (!a->rate_only) {
                            int donor = v >= 0 ? L : R;
                            double n = q->ion[s][donor], u = n > 0 ? q->momentum[s][donor] / n : 0;
                            F = v * n;
                            H = F * u;
                        }
                        if (axis == 2 && a->c.parallel && (!wall || open)) {
                            double nl = q->ion[s][L], nr = q->ion[s][R],
                                   ul = nl > 0 ? q->momentum[s][L] / nl : 0,
                                   ur = nr > 0 ? q->momentum[s][R] / nr : 0;
                            double rawl = nl, rawr = nr;
                            if (a->c.parallel_order == 2 && !wall) {
                                if (!a->rate_only) {
                                    nl = state(a, s, i, j, k - 1, 0, 1);
                                    nr = state(a, s, i, j, k, 0, -1);
                                }
                                ul = state(a, s, i, j, k - 1, 1, 1);
                                ur = state(a, s, i, j, k, 1, -1);
                            }
                            double mass = (s ? 32 : 30) * amu,
                                   cs = sqrt(kb * (a->c.chem_Ti + a->c.chem_Te) / mass),
                                   speed = fmax(fabs(ul), fabs(ur)) + cs;
                            double N = .5 * (nl * ul + nr * ur) - .5 * speed * (nr - nl), theta, cl,
                                   cr;
                            double ds = .5 * (length(a, i, j, k - 1) + length(a, i, j, k)),
                                   nu = .5 * (q->nu[s][L] + q->nu[s][R]);
                            double g = .5 * (a->g->cell[L].gravity_parallel +
                                             a->g->cell[R].gravity_parallel);
                            night_sg_coeff(a->c.chem_Ti, a->c.chem_Te, mass, nu, g, ds,
                                           q->electron[L], q->electron[R], &theta, &cl, &cr);
                            N = theta * N + cl * rawl - cr * rawr;
                            double area = f->parallel_area;
                            if (!a->rate_only) {
                                F += area * N;
                                H += area * N * (N >= 0 ? ul : ur);
                            }
                            rate += area * fmax(speed, theta * speed + fmax(cl, cr));
                        }
                        if (ownL) {
                            size_t t = L - off;
                            if (!a->rate_only) {
                                q->dn[s][t] -= F;
                                q->dm[s][t] -= H;
                            }
                            q->rates[s][t] += rate;
                        }
                        if (ownR) {
                            size_t t = R - off;
                            if (!a->rate_only) {
                                q->dn[s][t] += F;
                                q->dm[s][t] += H;
                            }
                            q->rates[s][t] += rate;
                        }
                        if (!a->rate_only && !left_exists && ownR)
                            boundary += F;
                        if (!a->rate_only && !right_exists && ownL)
                            boundary -= F;
                    }
                }
            }
        }
        q->boundary_rate += boundary;
    }
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) reduction(max : maximum_rate)
#endif
    for (int j = 2 * a->row_first; j < 2 * a->row_last; j++)
        for (int i = 0; i < a->fx; i++)
            for (int k = 0; k < a->ns; k++) {
                int g = id(a, i, j, k);
                size_t t = g - off;
                TubeCell *cell = a->g->cell + g;
                double pe = a->c.parallel && !a->rate_only
                                ? slope(a, q->electron, q->transport_log[2], i, j, k)
                                : 0;
                for (int s = 0; s < 2; s++) {
                    if (!a->rate_only) {
                        q->dn[s][t] /= cell->volume;
                        q->dm[s][t] /= cell->volume;
                        if (a->c.parallel) {
                            double pi = slope(a, q->ion[s], q->transport_log[s], i, j, k);
                            q->dm[s][t] +=
                                q->ion[s][g] * (cell->gravity_parallel -
                                                night_pressure_force(a->c.chem_Ti, a->c.chem_Te,
                                                                     (s ? 32 : 30) * amu, pi, pe));
                        } else
                            q->dm[s][t] = 0;
                    }
                    maximum_rate = fmax(maximum_rate, q->rates[s][t] / cell->volume);
                }
            }
    MPI_Allreduce(MPI_IN_PLACE, &maximum_rate, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &q->boundary_rate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    a->rate = fmax(a->rate, maximum_rate);
    q->transport_wall += wall3() - begin;
}
/**
 * Collectively verify finite/nonnegative molecular densities and consistent exact-vacuum
 * momenta.
 */
static int valid(A3 *a) {
    A3Night *q = a->night;
    size_t off = night_start(a), count = night_owned(a);
    int ok = 1;
    for (int s = 0; s < 2; s++)
        for (size_t t = off; t < off + count; t++)
            if (!isfinite(q->ion[s][t]) || q->ion[s][t] < 0 || !isfinite(q->momentum[s][t]) ||
                (q->ion[s][t] == 0 && q->momentum[s][t] != 0))
                ok = 0;
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return ok;
}
/**
 * Save the accepted pre-transport molecular state and first-stage boundary rate for the two-
 * stage update.
 */
void a3_night_transport_begin(A3 *a) {
    if (!a->night)
        return;
    A3Night *q = a->night;
    size_t off = night_start(a), bytes = night_owned(a) * sizeof(double);
    for (int s = 0; s < 2; s++) {
        memcpy(q->old_n[s], q->ion[s] + off, bytes);
        memcpy(q->old_m[s], q->momentum[s] + off, bytes);
    }
    q->boundary0 = q->boundary_rate;
}
/**
 * Apply the molecular first-stage density/ETD momentum predictor and cache its fixed-nu
 * coefficients for this exact attempted dt.
 */
int a3_night_predict(A3 *a, double dt) {
    if (!a->night)
        return 1;
    A3Night *q = a->night;
    size_t off = night_start(a), count = night_owned(a);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int s = 0; s < 2; s++)
        for (size_t t = 0; t < count; t++) {
            double decay, p1, p2;
            exponential_weights(q->nu[s][off + t], dt, &decay, &p1, &p2);
            q->etd_weight[s][0][t] = decay;
            q->etd_weight[s][1][t] = p1;
            q->etd_weight[s][2][t] = p2;
            q->first_m[s][t] = q->dm[s][t];
            q->ion[s][off + t] = q->old_n[s][t] + dt * q->dn[s][t];
            q->momentum[s][off + t] = decay * q->old_m[s][t] + dt * p1 * q->dm[s][t];
        }
    return valid(a);
}
/**
 * Complete the two-stage density and ETD2 momentum update using same-attempt coefficients, then
 * check positivity and stage boundary balance.
 */
int a3_night_correct(A3 *a, double dt) {
    if (!a->night)
        return 1;
    A3Night *q = a->night;
    size_t off = night_start(a), count = night_owned(a);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int s = 0; s < 2; s++)
        for (size_t t = 0; t < count; t++) {
            double decay = q->etd_weight[s][0][t], p1 = q->etd_weight[s][1][t],
                   p2 = q->etd_weight[s][2][t];
            q->ion[s][off + t] = .5 * q->old_n[s][t] + .5 * (q->ion[s][off + t] + dt * q->dn[s][t]);
            q->momentum[s][off + t] =
                decay * q->old_m[s][t] + dt * ((p1 - p2) * q->first_m[s][t] + p2 * q->dm[s][t]);
        }
    q->boundary_pending = .5 * dt * (q->boundary0 + q->boundary_rate);
    return valid(a);
}
/**
 * Commit the accepted molecular boundary-particle increment after the coupled step passes all
 * checks.
 */
void a3_night_accept(A3 *a) {
    if (a->night)
        a->night->boundary_total += a->night->boundary_pending;
}
