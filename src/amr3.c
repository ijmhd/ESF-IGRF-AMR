/*
 * AMR driver and conservative O+ fluxes. n is normalized by 1e12 m^-3; m=n*u_parallel. Private
 * levels -1/0/1 represent physical L0/L1/L2 in a three-level run.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only */
#include "amr3.h"
#include "night_transport.h"
#include "neutral_grid.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.14159265358979323846

/**
 * Return the nonnegative periodic index i modulo n.
 */
static int wrap(int i, int n) {
    i %= n;
    return i < 0 ? i + n : i;
}
/**
 * Map a canonical transverse cell to its backing leaf column; apply the configured
 * periodic/clamped indexing.
 */
static int cover(A3 *a, int i, int j) {
    if (a->c.periodic_x || a->c.amr_case == 1)
        i = wrap(i, a->fx);
    if (a->c.amr_case == 1)
        j = wrap(j, a->fy);
    if (i < 0 || i >= a->fx || j < 0 || j >= a->fy)
        return -1;
    return a->cover[j * a->fx + i];
}
/**
 * Test whether a backing column belongs to this rank, so shared faces update only local
 * conserved states.
 */
static int owned(A3 *a, int p) {
    return p >= a->owned_first && p < a->owned_first + a->owned_count;
}
/**
 * Read a transverse potential node using this solver's boundary/periodic indexing convention.
 */
static double ph(A3 *a, int i, int j) {
    return a->phi[j * (a->fx + 1) + i];
}
/**
 * Compute the two transverse potential derivatives at a geometric face and retain both metric
 * cross contributions.
 */
static void derivative(A3 *a, int axis, int i, int j, double d[2]) {
    int nx = a->fx, ny = a->fy;
    if (axis == 0) {
        int il = i ? i - 1 : i, ir = i < nx ? i + 1 : i, denom = ir - il;
        if (a->c.periodic_x && (i == 0 || i == nx)) {
            il = nx - 1;
            ir = 1;
            denom = 2;
            i = 0;
        }
        d[0] =
            nx * .5 * (ph(a, ir, j) + ph(a, ir, j + 1) - ph(a, il, j) - ph(a, il, j + 1)) / denom;
        d[1] = ny * (ph(a, i, j + 1) - ph(a, i, j));
    } else if (axis == 1) {
        int jl = j ? j - 1 : j, jr = j < ny ? j + 1 : j;
        d[0] = nx * (ph(a, i + 1, j) - ph(a, i, j));
        d[1] = ny * .5 * (ph(a, i, jr) + ph(a, i + 1, jr) - ph(a, i, jl) - ph(a, i + 1, jl)) /
               (jr - jl);
    } else {
        d[0] = nx * .5 * (ph(a, i + 1, j) + ph(a, i + 1, j + 1) - ph(a, i, j) - ph(a, i, j + 1));
        d[1] = ny * .5 * (ph(a, i, j + 1) + ph(a, i + 1, j + 1) - ph(a, i, j) - ph(a, i + 1, j));
    }
}
/**
 * Return the sign-consistent slope with smallest magnitude; return zero at an extremum.
 */
static double minmod(double a, double b, double c) {
    return a * b > 0 && a * c > 0 ? copysign(fmin(fabs(a), fmin(fabs(b), fabs(c))), a) : 0;
}
/**
 * Read the physical along-field segment length at the requested private refinement level.
 */
static double length_at(A3 *a, int level, int i, int j, int k) {
    k = k < 0 ? 0 : k >= a->ns ? a->ns - 1 : k;
    if (a->parallel_length[level])
        return a->parallel_length[level][(j * (a3_nx(a, level)) + i) * a->ns + k];
    int scale = a3_scale(level);
    double V = 0, A = 0;
    for (int v = 0; v < scale; v++)
        for (int u = 0; u < scale; u++) {
            int x = scale * i + u, y = scale * j + v;
            V += a->g->cell[tube_cell_id(a->g, x, y, k)].volume;
            A += .5 * (fabs(a->g->face[2][tube_face_id(a->g, 2, x, y, k)].parallel_area) +
                       fabs(a->g->face[2][tube_face_id(a->g, 2, x, y, k + 1)].parallel_area));
        }
    return V / A;
}
/**
 * Sample density or velocity with the caller's parallel boundary convention; used by the along-
 * field reconstruction.
 */
static double cell_value(A3 *a, int p, int k, int field) {
    int here = k < 0 ? 0 : k >= a->ns ? a->ns - 1 : k;
    A3Leaf l = a3_density_leaf(a, p, here);
    return a3_sample(a, field, l.level, l.i, l.j, k);
}
/**
 * Construct a bounded left/right parallel face state at an explicitly requested refinement
 * level.
 */
static double parallel_state_at(A3 *a, int level, int i, int j, int k, int field, int side) {
    double left = a3_sample(a, field, level, i, j, k - 1),
           center = a3_sample(a, field, level, i, j, k),
           right = a3_sample(a, field, level, i, j, k + 1);
    return center +
           .5 * side * minmod(1.3 * (center - left), .5 * (right - left), 1.3 * (right - center));
}
/**
 * Resolve the actual local density level, then reconstruct its parallel face state.
 */
static double parallel_state(A3 *a, int p, int k, int field, int side) {
    A3Leaf l = a3_density_leaf(a, p, k);
    return parallel_state_at(a, l.level, l.i, l.j, k, field, side);
}
/**
 * Resolve a backing cell's active density level before reading its physical parallel length.
 */
static double length(A3 *a, int p, int k) {
    A3Leaf l = a3_density_leaf(a, p, k);
    return length_at(a, l.level, l.i, l.j, k);
}
/**
 * Compute the nonuniform parallel logarithmic pressure gradient, using common Ne in the multi-
 * ion closure and exact-vacuum secants.
 */
static double log_pressure(A3 *a, int p, int k) {
    A3Leaf l = a3_density_leaf(a, p, k);
    int lo = k ? k - 1 : k, hi = k + 1 < a->ns ? k + 1 : k;
    double x = a3_sample(a, 0, l.level, l.i, l.j, k);
    if (x == 0)
        return 0;
    double nlo = a3_sample(a, 0, l.level, l.i, l.j, lo),
           nhi = a3_sample(a, 0, l.level, l.i, l.j, hi);
    double h = length_at(a, l.level, l.i, l.j, k),
           hl = .5 * (h + length_at(a, l.level, l.i, l.j, lo)),
           hr = .5 * (h + length_at(a, l.level, l.i, l.j, hi));
    double dl = nlo > 0 ? log(x / nlo) : (x - nlo) / x;
    double dr = nhi > 0 ? log(nhi / x) : (nhi - x) / x;
    if (a->night) {
        double ne = a3_night_electron(a, l.level, l.i, l.j, k),
               nel = a3_night_electron(a, l.level, l.i, l.j, lo),
               ner = a3_night_electron(a, l.level, l.i, l.j, hi);
        dl = (a->c.chem_Ti * dl + a->c.chem_Te * night_log_jump(nel, ne)) / a->c.temperature;
        dr = (a->c.chem_Ti * dr + a->c.chem_Te * night_log_jump(ne, ner)) / a->c.temperature;
    }
    if (k == 0)
        return dr / hr;
    if (k == a->ns - 1)
        return dl / hl;
    return (hr / hl * dl + hl / hr * dr) / (hl + hr);
}
/**
 * Evaluate the Bernoulli fitting function z/(exp(z)-1) with small-argument series and large-
 * argument limits.
 */
static double bernoulli(double z) {
    if (fabs(z) < 1e-4) {
        double z2 = z * z;
        return 1 - .5 * z + z2 / 12 - z2 * z2 / 720;
    }
    if (z > 50)
        return z * exp(-z) / (1 - exp(-z));
    if (z < -50)
        return -z / (1 - exp(z));
    return z / expm1(z);
}
/**
 * Prepare shared parallel flux/rate tables for a column-uniform density layout; segmented
 * layouts use a different path.
 */
static void parallel_cache(A3 *a, double cs) {
    if (!a->c.parallel)
        return;
    for (int slot = 0; slot < a->owned_count; slot++) {
        int p = a->global[slot];
        for (int k = 0; k <= a->ns; k++) {
            int t = slot * (a->ns + 1) + k, lo = k ? k - 1 : 0, hi = k < a->ns ? k : a->ns - 1;
            int wall = k == 0 || k == a->ns, open = a->c.open_ends;
            double raw_l = cell_value(a, p, lo, 0), raw_r = cell_value(a, p, hi, 0);
            double nl = raw_l, nr = raw_r, ml = cell_value(a, p, lo, 1),
                   mr = cell_value(a, p, hi, 1);
            double ul = nl > 0 ? ml / nl : 0, ur = nr > 0 ? mr / nr : 0;
            if (a->c.parallel_order == 2 && !wall) {
                nl = parallel_state(a, p, lo, 0, 1);
                nr = parallel_state(a, p, hi, 0, -1);
                if (a->c.parallel_relaxation) {
                    ul = parallel_state(a, p, lo, 2, 1);
                    ur = parallel_state(a, p, hi, 2, -1);
                    ml = nl * ul;
                    mr = nr * ur;
                } else {
                    ml = parallel_state(a, p, lo, 1, 1);
                    mr = parallel_state(a, p, hi, 1, -1);
                    ul = nl > 0 ? ml / nl : 0;
                    ur = nr > 0 ? mr / nr : 0;
                }
            }
            if (k == 0 && !open) {
                nl = nr;
                ml = -mr;
                ul = -ur;
            }
            if (k == a->ns && !open) {
                nr = nl;
                mr = -ml;
                ur = -ul;
            }
            double speed = fmax(fabs(ul), fabs(ur)) + cs;
            double N = .5 * (ml + mr) - .5 * speed * (nr - nl), M = 0, rate = speed;
            if (!a->c.parallel_relaxation)
                M = .5 * (ml * ul + mr * ur + cs * cs * (nl + nr)) - .5 * speed * (mr - ml);
            if (a->c.parallel_relaxation && (!wall || open)) {
                A3Metric *L = a->metric + p * a->ns + lo, *R = a->metric + p * a->ns + hi;
                double nu = .5 * (L->nu + R->nu), ds = .5 * (length(a, p, lo) + length(a, p, hi));
                if (nu > 0) {
                    double g = .5 * (L->gravity + R->gravity), kn = nu * ds / cs,
                           theta = 1 / (1 + kn * kn);
                    double pe = g * ds / (cs * cs), bl = bernoulli(-pe), br = bernoulli(pe);
                    double conductance = cs * cs / (nu * ds);
                    N = theta * N + (1 - theta) * conductance * (bl * raw_l - br * raw_r);
                    rate = fmax(rate, theta * rate + (1 - theta) * conductance * fmax(bl, br));
                }
                M = N * (N >= 0 ? ul : ur);
            }
            a->parallel_n[t] = wall && !open ? 0 : N;
            a->parallel_m[t] = M;
            a->parallel_rate[t] = rate;
        }
    }
}
/**
 * Locate the actual density owner at a canonical cell and along-field index; owner can differ
 * from electric-field backing.
 */
static int density_cover(A3 *a, int i, int j, int k) {
    int p = cover(a, i, j);
    if (p < 0)
        return -1;
    if (a->c.amr_case == 1)
        k = wrap(k, a->ns);
    else
        k = k < 0 ? 0 : k >= a->ns ? a->ns - 1 : k;
    return a3_owner(a, p, k);
}
/**
 * Append or merge a canonical face interval whose two active density owners and levels are
 * unchanged.
 */
static void append_face(A3 *a, int axis, int i, int j, int first, int last) {
    if (first >= last)
        return;
    int all = a->c.amr_case == 1, kl = first - (axis == 2), kr = first;
    if (axis == 2 && all) {
        kl = wrap(kl, a->ns);
        kr = wrap(kr, a->ns);
    }
    int left = density_cover(a, i - (axis == 0), j - (axis == 1), kl),
        right = density_cover(a, i, j, kr);
    if ((axis < 2 && left == right) || (!owned(a, left) && !owned(a, right)))
        return;
    int lk = kl < 0        ? 0
             : kl >= a->ns ? a->ns - 1
                           : kl,
        rk = kr < 0        ? 0
             : kr >= a->ns ? a->ns - 1
                           : kr;
    int ll = left >= 0 ? a3_density_leaf(a, left, lk).level : (a->c.amr_levels == 3 ? -1 : 0);
    int rl = right >= 0 ? a3_density_leaf(a, right, rk).level : (a->c.amr_levels == 3 ? -1 : 0);
    int level = ll > rl ? ll : rl, scale = a3_scale(level);
    if ((axis == 0 && j % scale) || (axis == 1 && i % scale) ||
        (axis == 2 && (i % scale || j % scale)))
        return;
    if (a->face_block_count == a->face_block_capacity) {
        a->face_block_capacity = a->face_block_capacity ? 2 * a->face_block_capacity : 1024;
        A3FaceBlock *next =
            realloc(a->face_blocks, (size_t)a->face_block_capacity * sizeof(A3FaceBlock));
        if (!next)
            a3_fail(a, "face topology allocation");
        a->face_blocks = next;
    }
    a->face_blocks[a->face_block_count++] =
        (A3FaceBlock){axis, i, j, first, last, left, right, level, ll, rl};
}
/**
 * Build conservative face blocks for the current segmented density mesh; identical interface
 * fluxes are shared by both adjacent volumes.
 */
void a3_face_topology(A3 *a) {
    free(a->face_blocks);
    a->face_blocks = NULL;
    a->face_block_count = a->face_block_capacity = 0;
    int all = a->c.amr_case == 1, slab = a->c.amr_slab ? a->c.amr_slab : a->ns;
    for (int axis = 0; axis < 3; axis++) {
        int ni = a->fx + (axis == 0 && !a->c.periodic_x && !all), nj = a->fy + (axis == 1 && !all);
        for (int row = 2 * a->row_first; row < 2 * a->row_last + (axis == 1); row++)
            for (int i = 0; i < ni; i++) {
                int j = row;
                if (j >= nj) {
                    if (axis == 1 && all && a->row_first > 0)
                        j = 0;
                    else
                        continue;
                }
                for (int k = 0; k < a->ns; k += slab) {
                    if (axis < 2)
                        append_face(a, axis, i, j, k, k + slab);
                    else {
                        append_face(a, axis, i, j, k, k + 1);
                        append_face(a, axis, i, j, k + 1, k + slab);
                    }
                }
                if (axis == 2 && !all)
                    append_face(a, axis, i, j, a->ns, a->ns + 1);
            }
    }
}
/**
 * Read or conservatively aggregate geometry for a requested physical density level, including
 * the macro level.
 */
static void metric_at(A3 *a, int level, int i, int j, int k, A3Metric *m) {
    memset(m, 0, sizeof(*m));
    int s = a3_scale(level);
    for (int v = 0; v < s; v++)
        for (int u = 0; u < s; u++) {
            TubeCell *g = a->g->cell + tube_cell_id(a->g, s * i + u, s * j + v, k);
            m->volume += g->volume;
            m->nu += g->volume * g->nu;
            m->gravity += g->volume * g->gravity_parallel;
        }
    m->nu /= m->volume;
    m->gravity /= m->volume;
}
/**
 * Cache static segment lengths and parallel relaxation factors for every supported level.
 */
void a3_parallel_geometry(A3 *a) {
    double cs = sqrt(1.380649e-23 * a->c.temperature / (16 * 1.66053906660e-27));
    for (int level = a->c.amr_levels == 3 ? -1 : 0; level < 2; level++) {
        int nx = a3_nx(a, level), ny = a3_ny(a, level);
        double *lengths = a3_alloc(a, (size_t)nx * ny * a->ns, sizeof(double));
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                for (int k = 0; k < a->ns; k++)
                    lengths[(j * nx + i) * a->ns + k] = length_at(a, level, i, j, k);
        a->parallel_length[level] = lengths;
        if (!a->c.amr_slab || !a->c.parallel_relaxation || a->c.night_chemistry)
            continue;
        double *response = a3_alloc(a, (size_t)nx * ny * (a->ns + 1) * 3, sizeof(double));
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                for (int k = 0; k <= a->ns; k++) {
                    int lo = k ? k - 1 : 0, hi = k < a->ns ? k : a->ns - 1;
                    A3Metric L, R;
                    metric_at(a, level, i, j, lo, &L);
                    metric_at(a, level, i, j, hi, &R);
                    double *r = response + 3 * ((j * nx + i) * (a->ns + 1) + k),
                           nu = .5 * (L.nu + R.nu);
                    r[0] = 1.;
                    if (nu > 0) {
                        double ds =
                            .5 * (length_at(a, level, i, j, lo) + length_at(a, level, i, j, hi));
                        double kn = nu * ds / cs, theta = 1 / (1 + kn * kn),
                               pe = .5 * (L.gravity + R.gravity) * ds / (cs * cs),
                               conductance = cs * cs / (nu * ds);
                        r[0] = theta;
                        r[1] = (1 - theta) * conductance * bernoulli(-pe);
                        r[2] = (1 - theta) * conductance * bernoulli(pe);
                    }
                }
        a->parallel_response[level] = response;
    }
}
/**
 * Evaluate parallel particle/momentum flux for one actual segmented interface; retain
 * collisional relaxation and electron-pressure fitting.
 */
static void segmented_parallel(A3 *a, int level, int i, int j, int k, double cs, double *N,
                               double *M, double *rate) {
    int lo = k ? k - 1 : 0, hi = k < a->ns ? k : a->ns - 1, wall = k == 0 || k == a->ns,
        open = a->c.open_ends;
    double raw_l = a3_sample(a, 0, level, i, j, lo), raw_r = a3_sample(a, 0, level, i, j, hi);
    double nl = raw_l, nr = raw_r, ml = a3_sample(a, 1, level, i, j, lo),
           mr = a3_sample(a, 1, level, i, j, hi);
    double ul = nl > 0 ? ml / nl : 0, ur = nr > 0 ? mr / nr : 0;
    if (a->c.parallel_order == 2 && !wall) {
        nl = parallel_state_at(a, level, i, j, lo, 0, 1);
        nr = parallel_state_at(a, level, i, j, hi, 0, -1);
        if (a->c.parallel_relaxation) {
            ul = parallel_state_at(a, level, i, j, lo, 2, 1);
            ur = parallel_state_at(a, level, i, j, hi, 2, -1);
            ml = nl * ul;
            mr = nr * ur;
        } else {
            ml = parallel_state_at(a, level, i, j, lo, 1, 1);
            mr = parallel_state_at(a, level, i, j, hi, 1, -1);
            ul = nl > 0 ? ml / nl : 0;
            ur = nr > 0 ? mr / nr : 0;
        }
    }
    if (k == 0 && !open) {
        nl = nr;
        ml = -mr;
        ul = -ur;
    }
    if (k == a->ns && !open) {
        nr = nl;
        mr = -ml;
        ur = -ul;
    }
    double speed = fmax(fabs(ul), fabs(ur)) + cs;
    *N = .5 * (ml + mr) - .5 * speed * (nr - nl);
    *M = 0;
    *rate = speed;
    if (!a->c.parallel_relaxation)
        *M = .5 * (ml * ul + mr * ur + cs * cs * (nl + nr)) - .5 * speed * (mr - ml);
    if (a->c.parallel_relaxation && (!wall || open)) {
        double mixed[3];
        const double *r;
        if (a->night) {
            A3Metric L, R;
            metric_at(a, level, i, j, lo, &L);
            metric_at(a, level, i, j, hi, &R);
            double ds = .5 * (length_at(a, level, i, j, lo) + length_at(a, level, i, j, hi));
            night_sg_coeff(a->c.chem_Ti, a->c.chem_Te, 16 * 1.66053906660e-27, .5 * (L.nu + R.nu),
                           .5 * (L.gravity + R.gravity), ds, a3_night_electron(a, level, i, j, lo),
                           a3_night_electron(a, level, i, j, hi), mixed, mixed + 1, mixed + 2);
            r = mixed;
        } else
            r = a->parallel_response[level] + 3 * ((j * (a3_nx(a, level)) + i) * (a->ns + 1) + k);
        *N = r[0] * (*N) + r[1] * raw_l - r[2] * raw_r;
        *rate = fmax(speed, r[0] * speed + fmax(r[1], r[2]));
        *M = (*N) * ((*N) >= 0 ? ul : ur);
    }
    if (wall && !open)
        *N = 0;
}
/**
 * Return a bounded material seed fraction from the advected marker and density; zero density
 * gives zero fraction.
 */
static double marker_fraction(A3 *a, int level, int axis, int i, int j, int k, int sign) {
    int q[3] = {i, j, k};
    q[axis] -= sign > 0;
    double n = a3_sample(a, 0, level, q[0], q[1], q[2]);
    return n > 0 ? fmax(0., fmin(1., a3_sample(a, 3, level, q[0], q[1], q[2]) / n)) : 0;
}
/**
 * Compute O+ finite-volume density/momentum RHS, shared face fluxes, gravity/pressure sources
 * and global stability/travel rates. Collective; sampling may populate lazy caches, so AMR face
 * loops remain on the main thread.
 */
void a3_transport(A3 *a, double *dn, double *dm) {
    double begin = wall3(), cs = sqrt(1.380649e-23 * a->c.temperature / (16 * 1.66053906660e-27));
    int n = a->owned_count * a->ns, periodic_all = a->c.amr_case == 1;
    memset(dn, 0, n * sizeof(double));
    memset(dm, 0, n * sizeof(double));
    if (a->marker_rhs)
        memset(a->marker_rhs, 0, n * sizeof(double));
    double *rate = a->rates;
    memset(rate, 0, n * sizeof(double));
    memset(a->transverse_rates, 0, n * sizeof(double));
    if (!a->c.amr_slab)
        parallel_cache(a, cs);
    if (a->audit_file)
        fputs(
            "axis,i,j,k,left_leaf,right_leaf,left_owned,right_owned,left_level,right_level,particle_flux,momentum_flux\n",
            a->audit_file);
    a->boundary_rate = 0;
    if (!periodic_all)
        for (int axis = 0; axis < 3; axis++) {
            int ni = a->fx + (axis == 0 && !a->c.periodic_x), nj = a->fy + (axis == 1);
            for (int j = 2 * a->row_first; j < 2 * a->row_last + (axis == 1); j++)
                if (j < nj)
                    for (int i = 0; i < ni; i++)
                        derivative(a, axis, i, j,
                                   a->potential_gradient[axis] + 2 * (j * (a->fx + 1) + i));
        }
    for (int face = 0; face < a->face_block_count; face++) {
        A3FaceBlock *f = a->face_blocks + face;
        int axis = f->axis, i = f->i, j = f->j, left = f->left, right = f->right, level = f->level,
            scale = a3_scale(level);
        int base[16], count = 0;
        double derivative_[16][2] = {{0}};
        for (int v = 0; v < (axis == 0 || axis == 2 ? scale : 1); v++)
            for (int u = 0; u < (axis == 1 || axis == 2 ? scale : 1); u++) {
                base[count] = tube_face_id(a->g, axis, i + u, j + v, 0);
                if (!periodic_all) {
                    const double *q =
                        a->potential_gradient[axis] + 2 * ((j + v) * (a->fx + 1) + i + u);
                    derivative_[count][0] = q[0];
                    derivative_[count][1] = q[1];
                }
                count++;
            }
        for (int k = f->first; k < f->last; k++) {
            int kl = k - (axis == 2), kr = k;
            if (axis == 2 && periodic_all) {
                kl = wrap(kl, a->ns);
                kr = wrap(kr, a->ns);
            }
            int wall = (axis == 0 && !a->c.periodic_x && !periodic_all && (i == 0 || i == a->fx)) ||
                       (axis == 1 && !periodic_all && (j == 0 || j == a->fy)) ||
                       (axis == 2 && !periodic_all && (k == 0 || k == a->ns));
            int open = (axis == 0 && a->c.regional_x) || (axis == 1 && a->c.open_apex) ||
                       (axis == 2 && a->c.open_ends);
            double positive = 0, negative = 0, area = 0;
            for (int q = 0; q < count; q++) {
                TubeFace *g = a->g->face[axis] + base[q] + k;
                double flux = periodic_all
                                  ? (axis == 0   ? .7 / (a->fy * a->ns)
                                     : axis == 1 ? .3 / (a->fx * a->ns)
                                                 : .2 / (a->fx * a->fy))
                                  : g->e[0] * derivative_[q][0] + g->e[1] * derivative_[q][1];
                if (wall && !open)
                    flux = 0;
                positive += fmax(flux, 0);
                negative += fmin(flux, 0);
                area += g->parallel_area;
            }
            double qabs = positive - negative, flux_n = 0, flux_m = 0, flux_p = 0,
                   parallel_rate = 0;
            for (int sign = -1; sign <= 1; sign += 2) {
                double q = sign > 0 ? positive : negative;
                if (q == 0)
                    continue;
                double density, second;
                a3_face_pair(a, level, axis, i / scale, j / scale, k, sign, &density, &second);
                flux_n += q * density;
                flux_m += a->c.parallel_relaxation ? q * density * second : q * second;
                if (a->marker)
                    flux_p += q * density *
                              marker_fraction(a, level, axis, i / scale, j / scale, k, sign);
            }
            if (axis == 2 && a->c.parallel) {
                double N, M, R;
                if (a->c.amr_slab)
                    segmented_parallel(a, level, i / scale, j / scale, k, cs, &N, &M, &R);
                else {
                    int t = a->slot[left] * (a->ns + 1) + k;
                    N = a->parallel_n[t];
                    M = a->parallel_m[t];
                    R = a->parallel_rate[t];
                }
                if (!(area > 0))
                    a3_fail(a, "nonpositive field-line face area");
                flux_n += area * N;
                flux_m += area * M;
                parallel_rate = area * R;
                if (a->marker && N != 0)
                    flux_p +=
                        area * N *
                        marker_fraction(a, level, axis, i / scale, j / scale, k, N > 0 ? 1 : -1);
            }
            if (owned(a, left) && kl >= 0 && kl < a->ns) {
                int t = a->slot[left] * a->ns + kl;
                dn[t] -= flux_n;
                dm[t] -= flux_m;
                rate[t] += qabs + parallel_rate;
                if (a->marker_rhs)
                    a->marker_rhs[t] -= flux_p;
                if (axis < 2)
                    a->transverse_rates[t] += qabs;
            }
            if (owned(a, right) && kr >= 0 && kr < a->ns) {
                int t = a->slot[right] * a->ns + kr;
                dn[t] += flux_n;
                dm[t] += flux_m;
                rate[t] += qabs + parallel_rate;
                if (a->marker_rhs)
                    a->marker_rhs[t] += flux_p;
                if (axis < 2)
                    a->transverse_rates[t] += qabs;
            }
            if (!periodic_all) {
                if (axis < 2) {
                    if (left < 0 && owned(a, right))
                        a->boundary_rate += flux_n;
                    if (right < 0 && owned(a, left))
                        a->boundary_rate -= flux_n;
                } else {
                    if (k == 0 && owned(a, right))
                        a->boundary_rate += flux_n;
                    if (k == a->ns && owned(a, left))
                        a->boundary_rate -= flux_n;
                }
            }
            if (a->audit_file && left >= 0 && right >= 0 && f->ll != f->rl)
                fprintf(a->audit_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.17g,%.17g\n", axis, i, j, k,
                        left, right, owned(a, left), owned(a, right), f->ll, f->rl, flux_n, flux_m);
        }
    }
    a->rate = 0;
    a->transverse_rate = 0;
    for (int p = 0; p < a->owned_count; p++)
        for (int k = 0; k < a->ns; k++) {
            int t = p * a->ns + k, id = a->global[p];
            A3Metric *m = a->metric + id * a->ns + k;
            if (!a3_active(a, id, k))
                continue;
            dn[t] /= m->volume;
            dm[t] /= m->volume;
            if (a->marker_rhs)
                a->marker_rhs[t] /= m->volume;
            if (a->c.parallel) {
                if (a->c.parallel_relaxation)
                    dm[t] += a->current_n[t] * (m->gravity - cs * cs * log_pressure(a, id, k));
                else
                    dm[t] += cs * cs * a->current_n[t] * (m->area[1] - m->area[0]) / m->volume +
                             a->current_n[t] * m->gravity;
            } else
                dm[t] = 0;
            a->rate = fmax(a->rate, rate[t] / m->volume);
            int level = a3_density_leaf(a, id, k).level;
            a->transverse_rate =
                fmax(a->transverse_rate, ldexp(1., -level) * a->transverse_rates[t] / m->volume);
        }
    MPI_Allreduce(MPI_IN_PLACE, &a->rate, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &a->boundary_rate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &a->transverse_rate, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    a->t_transport += wall3() - begin;
}

/**
 * Execute the complete current-stage chain: prepare/halo, integrated potential solve, O+
 * transport, then molecular transport. The potential always belongs to this stage's density and
 * chemistry state.
 */
static void rhs(A3 *a, double *n, double *m, double *dn, double *dm, double time) {
    a->rhs_calls++;
    a->time = time;
    a3_prepare(a, n, m);
    if (a->c.amr_case != 1)
        a3_potential(a);
    a3_transport(a, dn, dm);
    a3_night_transport(a);
}
/**
 * Flatten node coordinates with the along-field coordinate contiguous.
 */
static int node_id(TubeGrid *g, int i, int j, int k) {
    return (j * (g->c.nx + 1) + i) * (g->c.ns + 1) + k;
}
/**
 * Build Cartesian synthetic geometry for transport/potential verification rather than tracing a
 * physical magnetic field.
 */
static void make_box(TubeGrid *g, GridConfig c) {
    memset(g, 0, sizeof(*g));
    g->c = c;
    int nx = c.nx, ny = c.ny, ns = c.ns;
    g->node = calloc((size_t)(nx + 1) * (ny + 1) * (ns + 1), sizeof(V3));
    g->cell = calloc((size_t)nx * ny * ns, sizeof(TubeCell));
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++)
            for (int k = 0; k <= ns; k++)
                g->node[node_id(g, i, j, k)] = (V3){(double)i / nx, (double)j / ny, (double)k / ns};
    g->min_volume = 1. / (nx * ny * ns);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            for (int k = 0; k < ns; k++) {
                TubeCell *p = g->cell + tube_cell_id(g, i, j, k);
                p->volume = g->min_volume;
                p->height = 438e3;
                p->B = 1;
                p->b = (V3){0, 0, 1};
                p->gx = (V3){1, 0, 0};
                p->gy = (V3){0, 1, 0};
                p->pF = 1;
                p->position = (V3){(i + .5) / nx, (j + .5) / ny, (k + .5) / ns};
            }
    for (int a = 0; a < 3; a++) {
        int count = tube_face_count(g, a);
        g->face[a] = calloc(count, sizeof(TubeFace));
        for (int p = 0; p < count; p++) {
            if (a == 0)
                g->face[a][p].e[1] = -1. / (ny * ns);
            if (a == 1)
                g->face[a][p].e[0] = 1. / (nx * ns);
            if (a == 2)
                g->face[a][p].parallel_area = 1. / (nx * ny);
        }
    }
}
/**
 * Collect disjoint computed geometry row ranges into replicated canonical geometry, in MPI-
 * count-safe chunks.
 */
static void gather_geometry(A3 *a, void *data, int rows, size_t rowbytes, const int *cuts) {
    if ((size_t)rows * rowbytes > 2147483647)
        a3_fail(a, "geometry collective byte count overflow");
    int *counts = a3_alloc(a, 2 * a->size, sizeof(int)), *offset = counts + a->size;
    for (int r = 0; r < a->size; r++) {
        int lo = cuts ? cuts[r] : r * rows / a->size,
            hi = cuts ? cuts[r + 1] : (r + 1) * rows / a->size;
        counts[r] = (int)((hi - lo) * rowbytes);
        offset[r] = (int)(lo * rowbytes);
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_BYTE, data, counts, offset, MPI_BYTE, MPI_COMM_WORLD);
    free(counts);
}
/**
 * Distribute tube tracing and cell/face metric construction among ranks, then assemble the
 * replicated geometry and diagnostics.
 */
static void parallel_geometry(A3 *a, TubeGrid *g, Igrf *field, GridConfig fine) {
    if (tube_geometry_allocate(g, fine))
        a3_fail(a, "canonical geometry allocation");
    int *cuts = a3_alloc(a, a->size + 1, sizeof(int));
    double total = 0, prefix = 0;
    for (int j = 0; j <= fine.ny; j++)
        total += sqrt(fine.hmin + (fine.hmax - fine.hmin) * j / fine.ny);
    int next = 1;
    for (int j = 0; j <= fine.ny && next < a->size; j++) {
        prefix += sqrt(fine.hmin + (fine.hmax - fine.hmin) * j / fine.ny);
        if (prefix >= next * total / a->size) {
            cuts[next++] = j + 1;
        }
    }
    cuts[a->size] = fine.ny + 1;
    double begin = wall3();
    for (int j = cuts[a->rank]; j < cuts[a->rank + 1]; j++)
        for (int i = 0; i <= fine.nx; i++) {
            V3 *line = g->node + ((size_t)j * (fine.nx + 1) + i) * (fine.ns + 1);
            if (tube_trace(field, fine.lon + fine.span * ((double)i / fine.nx - .5),
                           fine.hmin + (fine.hmax - fine.hmin) * j / fine.ny, fine.hbase,
                           fine.trace_step, fine.ns, line))
                a3_fail(a, "parallel field-line tracing");
        }
    fprintf(stdout, "Geometry MPI rank %d traced %d node rows in %.3fs\n", a->rank,
            cuts[a->rank + 1] - cuts[a->rank], wall3() - begin);
    fflush(stdout);
    gather_geometry(a, g->node, fine.ny + 1, (size_t)(fine.nx + 1) * (fine.ns + 1) * sizeof(V3),
                    cuts);
    free(cuts);
    int first = a->rank * fine.ny / a->size, last = (a->rank + 1) * fine.ny / a->size;
    if (tube_geometry_rows(g, field, first, last))
        a3_fail(a, "parallel canonical cell/face metrics");
    gather_geometry(a, g->cell, fine.ny, (size_t)fine.nx * fine.ns * sizeof(TubeCell), NULL);
    for (int axis = 0; axis < 3; axis++) {
        cuts = a3_alloc(a, a->size + 1, sizeof(int));
        for (int r = 0; r <= a->size; r++)
            cuts[r] = r * fine.ny / a->size + (axis == 1 && r == a->size);
        gather_geometry(
            a, g->face[axis], fine.ny + (axis == 1),
            (size_t)(fine.nx + (axis == 0)) * (fine.ns + (axis == 2)) * sizeof(TubeFace), cuts);
        free(cuts);
    }
    double stats[] = {g->max_alignment_error, g->max_parallel_electric_metric_error,
                      g->max_magnetic_flux_jump};
    MPI_Allreduce(MPI_IN_PLACE, stats, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &g->min_volume, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    g->max_alignment_error = stats[0];
    g->max_parallel_electric_metric_error = stats[1];
    g->max_magnetic_flux_jump = stats[2];
}
/**
 * Load and validate a geometry cache, or construct geometry and optionally write a new cache;
 * initialize conductivity responses afterward.
 */
static void geometry(A3 *status, TubeGrid *g, Config3 c) {
    GridConfig fine = c.grid;
    fine.nx *= 2;
    fine.ny *= 2;
    if (c.amr_case == 1 || c.amr_case == 2) {
        make_box(g, fine);
        return;
    }
    Igrf field;
    if (solver3_load_field(&c, &field))
        a3_fail(status, "invalid magnetic field");
    int cached = 0;
    if (!status->rank)
        cached = a3_geometry_cache(status, g, c, 0);
    MPI_Bcast(&cached, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!cached) {
        if (!status->rank) {
            fprintf(stdout, "Building canonical %dx%dx%d geometry using %d MPI ranks\n", fine.nx,
                    fine.ny, fine.ns, status->size);
            fflush(stdout);
        }
        parallel_geometry(status, g, &field, fine);
        if (!status->rank) {
            tube_conductivity(g, c.e_peak);
            if (c.collision_profile[0] && tube_conductivity_table(g, c.e_peak, c.collision_profile))
                a3_fail(status, "invalid AMR collision table");
            a3_geometry_cache(status, g, c, 1);
        }
        size_t bytes = (size_t)fine.nx * fine.ny * fine.ns * sizeof(TubeCell);
        if (bytes > 2147483647)
            a3_fail(status, "cell geometry byte overflow");
        MPI_Bcast(g->cell, (int)bytes, MPI_BYTE, 0, MPI_COMM_WORLD);
        return;
    }
    if (!status->rank) {
    } else {
        memset(g, 0, sizeof(*g));
        g->c = fine;
        g->node =
            a3_alloc(status, (size_t)(fine.nx + 1) * (fine.ny + 1) * (fine.ns + 1), sizeof(V3));
        g->cell = a3_alloc(status, (size_t)fine.nx * fine.ny * fine.ns, sizeof(TubeCell));
        for (int axis = 0; axis < 3; axis++)
            g->face[axis] = a3_alloc(status, tube_face_count(g, axis), sizeof(TubeFace));
    }
    size_t bytes[] = {((size_t)fine.nx + 1) * (fine.ny + 1) * (fine.ns + 1) * sizeof(V3),
                      (size_t)fine.nx * fine.ny * fine.ns * sizeof(TubeCell)};
    if (bytes[0] > 2147483647 || bytes[1] > 2147483647)
        a3_fail(status, "geometry exceeds MPI byte count");
    MPI_Bcast(g->node, (int)bytes[0], MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g->cell, (int)bytes[1], MPI_BYTE, 0, MPI_COMM_WORLD);
    for (int axis = 0; axis < 3; axis++) {
        size_t n = (size_t)tube_face_count(g, axis) * sizeof(TubeFace);
        if (n > 2147483647)
            a3_fail(status, "face geometry exceeds MPI byte count");
        MPI_Bcast(g->face[axis], (int)n, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
}
/**
 * Open a named output in the current run directory; failures are collective-fatal through
 * a3_fail.
 */
static FILE *output(A3 *a, const char *name, const char *mode) {
    char path[1500];
    snprintf(path, sizeof(path), "%s/%s", a->c.out, name);
    FILE *f = fopen(path, mode);
    if (!f)
        a3_fail(a, "cannot write AMR output");
    return f;
}
/**
 * Write effective geometry, numerical controls, units and binary layout metadata for this AMR
 * run.
 */
static void metadata(A3 *a) {
    FILE *f = output(a, "metadata.json", "w");
    fprintf(
        f,
        "{\n\"model\":\"field-aligned ESF with composite anisotropic AMR\",\"amr_mode\":%d,\"case\":%d,\n",
        a->c.amr_mode, a->c.amr_case);
    fprintf(
        f,
        "\"amr_physical_levels\":%d,\"logical_root_nx\":%d,\"logical_root_ny\":%d,\"private_level_offset\":%d,\n",
        a->c.amr_levels, a->c.amr_levels == 3 ? a->nx / 2 : a->nx,
        a->c.amr_levels == 3 ? a->ny / 2 : a->ny, a->c.amr_levels == 3 ? 1 : 0);
    fprintf(f, "\"seed_apex_longitude_deg\":%.17g,\"seed_longitude_offset_deg\":%.17g,\n",
            a->c.grid.lon + a->c.seed_offset_deg, a->c.seed_offset_deg);
    fprintf(f, "\"night_chemistry\":%d,\"night_reactions\":%d,\"utc_start_unix_s\":%.17g,\n",
            a->c.night_chemistry, a->c.night_reactions, a->c.utc_start);
    fprintf(f, "\"root_nx\":%d,\"root_ny\":%d,\"fine_nx\":%d,\"fine_ny\":%d,\"ns\":%d,\n", a->nx,
            a->ny, a->fx, a->fy, a->ns);
    fprintf(f, "\"leaf_columns\":%d,\"refined_roots\":%d,\"potential_unknowns\":%d,\"ranks\":%d,\n",
            a->nleaf, a->split, a->ndof, a->size);
    fprintf(f, "\"field\":\"%s\",\"year\":%.17g,\"dipole_equator_nT\":%.17g,\n",
            a->c.aligned_dipole ? "aligned-dipole"
            : a->c.degree == 1  ? "tilted-dipole"
                                : "IGRF-14",
            a->c.year, a->c.dipole_equator_nT);
    fprintf(
        f,
        "\"periodic_x\":%d,\"E_peak_m3\":%.17g,\"exb_scheme\":\"%s\",\"beta\":%.17g,\"cfl\":%.17g,\n",
        a->c.periodic_x, a->c.e_peak, a->c.exb_pdm4 ? "PDM4 volume moments" : "MUSCL",
        a->c.pdm_beta, a->c.cfl);
    fprintf(
        f,
        "\"temperature_K\":%.17g,\"seed_amplitude\":%.17g,\"seed_height_m\":%.17g,\"seed_width_m\":%.17g,\n",
        a->c.temperature, a->c.amplitude, a->c.seed_height, a->c.seed_width);
    fprintf(
        f,
        "\"apex_min_m\":%.17g,\"apex_max_m\":%.17g,\"longitude_center_deg\":%.17g,\"longitude_span_deg\":%.17g,\n",
        a->c.grid.hmin, a->c.grid.hmax, a->c.grid.lon, a->c.grid.span);
    fprintf(
        f,
        "\"foot_height_m\":%.17g,\"trace_step_m\":%.17g,\"parallel\":%d,\"parallel_order\":%d,\"parallel_relaxation\":%d,\"hall\":%d,\"rtol\":%.17g,\"dt_max_s\":%.17g,\"background_profile\":%d,\n",
        a->c.grid.hbase, a->c.grid.trace_step, a->c.parallel, a->c.parallel_order,
        a->c.parallel_relaxation, a->c.hall, a->c.rtol, a->c.dtmax, a->c.background_profile);
    fprintf(
        f,
        "\"regional_x\":%d,\"open_apex\":%d,\"open_ends\":%d,\"seed_east_width_m\":%.17g,\"collision_table\":%s,\"dynamic_regrid\":%s,\n",
        a->c.regional_x, a->c.open_apex, a->c.open_ends, a->c.seed_east,
        a->c.collision_profile[0] ? "true" : "false", a->c.amr_mode == 3 ? "true" : "false");
    fprintf(f, "\"cached_COO_assembly\":%s,\"partition_canonical_work_weight\":%.17g,\n",
            a->c.amr_coo ? "true" : "false", a->c.amr_work_weight);
    fprintf(
        f,
        "\"amr_interval_s\":%.17g,\"amr_buffer_roots\":%d,\"amr_hold_checks\":%d,\"depletion_ratio_on\":%.17g,\"depletion_ratio_off\":%.17g,\"gradient_on\":%.17g,\"gradient_off\":%.17g,\n",
        a->c.amr_interval, a->c.amr_buffer, a->c.amr_hold, a->c.amr_deplete_on,
        a->c.amr_deplete_off, a->c.amr_grad_on, a->c.amr_grad_off);
    fputs(
        "\"integrator\":\"ETD2, SSPRK2 density, common timestep\",\"subcycling\":false,\n"
        "\"potential\":\"constrained conforming Q1; complete field-line Pedersen/Hall tensor; PETSc FGMRES/GAMG with true residual\",\n"
        "\"state_format\":\"little-endian float64 [2,leaf_columns,ns], density and parallel momentum\",\n"
        "\"potential_format\":\"float64 [fine_ny+1,fine_nx+1] constrained Q1 potential\",\n"
        "\"canonical_geometry_format\":\"float64 [fine_ny,fine_nx,ns,6], x,y,z,height,B,volume\",\n"
        "\"canonical_nodes_format\":\"float64 [fine_ny+1,fine_nx+1,ns+1,3], ECEF node positions\",\n"
        "\"leaf_metrics_format\":\"float64 [leaf_columns,ns,6], volume,bg,nu,gravity,area_lo,area_hi\",\n"
        "\"geometry\":\"replicated canonical fine subcells; coarse volumes and shared faces are exact sums\",\n"
        "\"mpi\":\"evolved state in owned+ghost slots; three-root-row halos; PETSc distributed rows; state gathered only for output; periodic tracer uses full peer slabs\",\n"
        "\"physical_boundaries\":\"zero polarization at apex edges; regional longitude d_xi(phi)=0 with open particles even at zero imposed E; configured open ends/apex; periodic longitude when selected\",\n"
        "\"physical_boundary_reconstruction\":\"MUSCL in two physical layers; periodic stencil elsewhere\",\n"
        "\"momentum_prolongation\":\"limited velocity centered by child particle weights, then conservative density times velocity\",\n"
        "\"timer\":\"macOS uptime excludes system sleep; monotonic clock on other platforms\"\n}\n",
        f);
    fclose(f);
}
/**
 * Collect owned states, reconstruct canonical display density, and write time/mesh-bound outputs
 * plus the native checkpoint. Display reconstruction is not an evolved finest-grid solution.
 */
static void snapshot(A3 *a, double *n, double *m, int frame, double time, int step, double mass0,
                     const int *quiet, double next_regrid, double travel, int retries) {
    int *counts = a3_alloc(a, a->size, sizeof(int)), *offsets = a3_alloc(a, a->size, sizeof(int));
    for (int r = 0; r < a->size; r++) {
        offsets[r] = a->row[a->partition[r]] * a->ns;
        counts[r] = (a->row[a->partition[r + 1]] - a->row[a->partition[r]]) * a->ns;
    }
    double mass = a3_mass(a, n), minimum = INFINITY, maximum = -INFINITY;
    for (int p = 0; p < a->owned_count; p++)
        for (int k = 0; k < a->ns; k++)
            if (a3_active(a, a->global[p], k)) {
                minimum = fmin(minimum, n[p * a->ns + k]);
                maximum = fmax(maximum, n[p * a->ns + k]);
            }
    MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    double *N = NULL, *M = NULL, *P = NULL;
    if (!a->rank && a->marker)
        P = a3_alloc(a, (size_t)a->nleaf * a->ns, sizeof(double));
    if (!a->rank) {
        N = a3_alloc(a, (size_t)a->nleaf * a->ns, sizeof(double));
        M = a3_alloc(a, (size_t)a->nleaf * a->ns, sizeof(double));
    }
    MPI_Gatherv(n, a->owned_count * a->ns, MPI_DOUBLE, N, counts, offsets, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    MPI_Gatherv(m, a->owned_count * a->ns, MPI_DOUBLE, M, counts, offsets, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    if (a->marker)
        MPI_Gatherv(a->marker, a->owned_count * a->ns, MPI_DOUBLE, P, counts, offsets, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
    if (a->c.amr_structure) {
        int count = 2 * (a->row_last - a->row_first) * a->fx * a->ns;
        double *display = a3_alloc(a, count, sizeof(double)), *all = NULL;
        for (int j = 2 * a->row_first; j < 2 * a->row_last; j++)
            for (int i = 0; i < a->fx; i++)
                for (int k = 0; k < a->ns; k++)
                    display[((j - 2 * a->row_first) * a->fx + i) * a->ns + k] =
                        a3_sample(a, 0, 1, i, j, k);
        int *c = a3_alloc(a, 2 * a->size, sizeof(int)), *o = c + a->size;
        for (int r = 0; r < a->size; r++) {
            o[r] = 2 * a->partition[r] * a->fx * a->ns;
            c[r] = 2 * (a->partition[r + 1] - a->partition[r]) * a->fx * a->ns;
        }
        if (!a->rank)
            all = a3_alloc(a, (size_t)a->fx * a->fy * a->ns, sizeof(double));
        MPI_Gatherv(display, count, MPI_DOUBLE, all, c, o, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (!a->rank) {
            char name[80];
            snprintf(name, sizeof(name), "display_density_%04d.bin", frame);
            FILE *f = output(a, name, "wb");
            size_t size = (size_t)a->fx * a->fy * a->ns;
            if (fwrite(all, sizeof(double), size, f) != size)
                a3_fail(a, "display reconstruction write");
            fclose(f);
        }
        free(display);
        free(all);
        free(c);
    }
    free(counts);
    free(offsets);
    if (a->rank)
        return;
    if (!frame) {
        metadata(a);
        if (a->c.collision_profile[0]) {
            FILE *src = fopen(a->c.collision_profile, "rb"),
                 *dst = output(a, "collision_profile.csv", "wb");
            if (!src)
                a3_fail(a, "collision input copy");
            unsigned char block[4096];
            size_t size;
            while ((size = fread(block, 1, sizeof(block), src)))
                if (fwrite(block, 1, size, dst) != size)
                    a3_fail(a, "collision copy write");
            if (ferror(src))
                a3_fail(a, "collision copy read");
            fclose(src);
            fclose(dst);
        }
        FILE *f = output(a, "leaves.csv", "w");
        fputs("id,level,i,j,root\n", f);
        for (int p = 0; p < a->nleaf; p++) {
            A3Leaf l = a->leaf[p];
            fprintf(f, "%d,%d,%d,%d,%d\n", p, l.level, l.i, l.j, l.root);
        }
        fclose(f);
        f = output(a, "geometry_cells.bin", "wb");
        for (int id = 0; id < a->fx * a->fy * a->ns; id++) {
            TubeCell *c = a->g->cell + id;
            double values[] = {c->position.x, c->position.y, c->position.z,
                               c->height,     c->B,          c->volume};
            if (fwrite(values, sizeof(double), 6, f) != 6)
                a3_fail(a, "canonical geometry write");
        }
        fclose(f);
        f = output(a, "nodes.bin", "wb");
        for (int id = 0; id < (a->fx + 1) * (a->fy + 1) * (a->ns + 1); id++) {
            V3 c = a->g->node[id];
            double values[] = {c.x, c.y, c.z};
            if (fwrite(values, sizeof(double), 3, f) != 3)
                a3_fail(a, "canonical node write");
        }
        fclose(f);
        f = output(a, "leaf_metrics.bin", "wb");
        for (int p = 0; p < a->nleaf * a->ns; p++) {
            A3Metric c = a->metric[p];
            double values[] = {c.volume, c.bg, c.nu, c.gravity, c.area[0], c.area[1]};
            if (fwrite(values, sizeof(double), 6, f) != 6)
                a3_fail(a, "metric write");
        }
        fclose(f);
        f = output(a, "diagnostics.csv", "w");
        fputs(
            "step,time_s,min_density,max_density,mass,relative_mass_balance,potential_residual,cfl_rate,solves,wall_s,integrated_boundary_particle_flux\n",
            f);
        fclose(f);
        f = output(a, "frames.csv", "w");
        fputs("frame,time_s,state,potential,mesh,leaf_metrics\n", f);
        fclose(f);
    }
    if (!frame || a->output_mesh_epoch != a->mesh_epoch) {
        char mesh[80];
        snprintf(mesh, sizeof(mesh), "leaves_%04d.csv", a->mesh_epoch);
        FILE *f = output(a, mesh, "w");
        fputs("id,level,i,j,root\n", f);
        for (int p = 0; p < a->nleaf; p++) {
            A3Leaf l = a->leaf[p];
            fprintf(f, "%d,%d,%d,%d,%d\n", p, l.level, l.i, l.j, l.root);
        }
        fclose(f);
        snprintf(mesh, sizeof(mesh), "leaf_metrics_%04d.bin", a->mesh_epoch);
        f = output(a, mesh, "wb");
        for (int p = 0; p < a->nleaf * a->ns; p++) {
            A3Metric c = a->metric[p];
            double values[] = {c.volume, c.bg, c.nu, c.gravity, c.area[0], c.area[1]};
            if (fwrite(values, sizeof(double), 6, f) != 6)
                a3_fail(a, "mesh metric write");
        }
        fclose(f);
        a->output_mesh_epoch = a->mesh_epoch;
    }
    char name[120];
    snprintf(name, sizeof(name), "state_%04d.bin", frame);
    FILE *f = output(a, name, "wb");
    size_t count = (size_t)a->nleaf * a->ns;
    if (fwrite(N, sizeof(double), count, f) != count ||
        fwrite(M, sizeof(double), count, f) != count)
        a3_fail(a, "state write");
    fclose(f);
    if (a->phi) {
        snprintf(name, sizeof(name), "phi_%04d.bin", frame);
        f = output(a, name, "wb");
        if (fwrite(a->phi, sizeof(double), (a->fx + 1) * (a->fy + 1), f) !=
            (size_t)(a->fx + 1) * (a->fy + 1))
            a3_fail(a, "potential write");
        fclose(f);
    }
    f = output(a, "frames.csv", "a");
    fprintf(f, "%d,%.17g,state_%04d.bin,", frame, time, frame);
    if (a->phi)
        fprintf(f, "phi_%04d.bin", frame);
    fprintf(f, ",leaves_%04d.csv,leaf_metrics_%04d.bin\n", a->mesh_epoch, a->mesh_epoch);
    fclose(f);
    f = output(a, "diagnostics.csv", "a");
    fprintf(f, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.6f,%.17g\n", step, time, minimum,
            maximum, mass, (mass - mass0 - a->boundary_total + a3_night_loss(a)) / mass0,
            a->potential_residual, a->rate, a->solves, wall3() - a->wall0, a->boundary_total);
    fclose(f);
    printf(
        "AMR t=%.3f step=%d columns=%d cells=%d potential=%d min=%.3g residual=%.2g mass=%.2g wall=%.1fs\n",
        time, step, a->nleaf, a->nleaf * a->ns, a->ndof, minimum, a->potential_residual,
        (mass - mass0 - a->boundary_total + a3_night_loss(a)) / mass0, wall3() - a->wall0);
    fflush(stdout);
    a3_structure_output(a, N, P, frame, time);
    a3_macro_output(a, frame, time);
    a3_checkpoint_write(a, N, M, P, time, step, mass0, quiet, next_regrid, travel, retries);
    free(N);
    free(M);
    free(P);
}
/**
 * Exercise constant/background/smooth conservative transfers and validate particle/momentum
 * invariants in the dedicated remap case.
 */
static void remap_tests(A3 *base, TubeGrid *g) {
    double worst = 0, field_error = 0, velocity_error = 0;
    for (int test = 1; test <= 4; test++) {
        A3 old;
        a3_mesh(&old, g, base->c, base->rank, base->size, 0);
        double *n = a3_alloc(&old, (size_t)old.slots * old.ns, sizeof(double)),
               *m = a3_alloc(&old, (size_t)old.slots * old.ns, sizeof(double));
        a3_initialize(&old, n, m, test);
        double mass = a3_mass(&old, n), momentum = a3_mass(&old, m), momentum_norm = 0;
        for (int p = 0; p < old.owned_count; p++)
            for (int k = 0; k < old.ns; k++)
                momentum_norm +=
                    fabs(m[p * old.ns + k]) * old.metric[old.global[p] * old.ns + k].volume;
        MPI_Allreduce(MPI_IN_PLACE, &momentum_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int iteration = 0; iteration < 6; iteration++) {
            A3 next;
            a3_mesh(&next, g, base->c, base->rank, base->size, (iteration + 1) % 3);
            double *N = a3_alloc(&next, (size_t)next.slots * next.ns, sizeof(double)),
                   *M = a3_alloc(&next, (size_t)next.slots * next.ns, sizeof(double));
            a3_remap(&old, &next, n, m, N, M);
            worst = fmax(worst, fabs(a3_mass(&next, N) - mass) / mass);
            worst = fmax(worst, fabs(a3_mass(&next, M) - momentum) / fmax(momentum_norm, 1e-300));
            for (int p = 0; p < next.owned_count; p++)
                for (int k = 0; k < next.ns; k++) {
                    int id = p * next.ns + k;
                    if (!a3_active(&next, next.global[p], k))
                        continue;
                    if (test < 4)
                        velocity_error = fmax(velocity_error, fabs(M[id] / N[id] - 2.5));
                    if (test <= 2) {
                        double truth =
                            test == 1 ? 1. : next.metric[next.global[p] * next.ns + k].bg;
                        double err = fabs(N[id] / truth - 1);
                        if (err > field_error && err > 1e-10 && next.c.amr_levels == 3)
                            fprintf(
                                stderr,
                                "remap equilibrium diagnostic test=%d cycle=%d p=%d k=%d level=%d N=%.17g truth=%.17g error=%.17g\n",
                                test, iteration, next.global[p], k,
                                a3_density_leaf(&next, next.global[p], k).level, N[id], truth, err);
                        if (err > field_error && err > 1e-10 && next.c.amr_levels == 3) {
                            A3Leaf l = a3_density_leaf(&next, next.global[p], k);
                            int fi = a3_scale(l.level) * l.i, fj = a3_scale(l.level) * l.j,
                                b = (fj / 2) * old.nx + fi / 2, macro = a3_macro_index(&old, b),
                                t = b * old.ns + k;
                            fprintf(
                                stderr,
                                "DETAIL old_macro N %.17g bg %.17g ratio %.17g L1 N %.17g bg %.17g A %.17g sx %.17g sy %.17g ready %d\n",
                                old.macro_state[0][macro * old.ns + k],
                                old.macro_bg[macro * old.ns + k],
                                old.macro_state[0][macro * old.ns + k] /
                                    old.macro_bg[macro * old.ns + k],
                                old.base_n[t], old.root_bg[t], old.amp_n[t], old.sx_n[t],
                                old.sy_n[t], old.ready[t]);
                        }
                        field_error = fmax(field_error, err);
                    }
                }
            free(n);
            free(m);
            a3_free(&old);
            old = next;
            a3_rebase(&old);
            n = N;
            m = M;
        }
        free(n);
        free(m);
        a3_free(&old);
    }
    MPI_Allreduce(MPI_IN_PLACE, &worst, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &field_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &velocity_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (worst > 1e-11 || field_error > 1e-10 || velocity_error > 1e-10) {
        fprintf(stderr, "remap errors mass=%.17g equilibrium=%.17g velocity=%.17g\n", worst,
                field_error, velocity_error);
        a3_fail(base, "physical-volume transfer invariant failed");
    }
    if (!base->rank) {
        FILE *f = output(base, "remap_checks.json", "w");
        fprintf(
            f,
            "{\"cycles_per_field\":6,\"fields\":\"constant density, Chapman equilibrium, deep depletion, signed momentum\",\"max_mass_and_momentum_relative_error\":%.17g,\"max_constant_or_equilibrium_relative_error\":%.17g,\"max_constant_velocity_error\":%.17g}\n",
            worst, field_error, velocity_error);
        fclose(f);
        printf(
            "PASS 3D volume remaps: mass/momentum %.3g, equilibrium %.3g, constant velocity %.3g\n",
            worst, field_error, velocity_error);
    }
}
/**
 * Import a uniform-grid restart by preserving original coarse-cell particle and signed momentum
 * integrals on the canonical AMR geometry.
 */
static void import_static(A3 *a, TubeGrid *g, double *n, double *m, double *time, int *step,
                          double *mass0) {
    Config3 cfg = a->c;
    cfg.amr_mode = -1;
    cfg.exb_pdm4 = 0;
    cfg.petsc_potential = 0;
    /* Every original coarse node is already present in the canonical geometry.
     * Rebuild original coarse trilinear volumes from those nodes, rather than
     * tracing the same field lines again or allocating an unused static solver. */
    Solver3 src = {.c = cfg, .rank = a->rank, .size = a->size, .wall0 = wall3()};
    if (solver3_load_field(&cfg, &src.field) || tube_geometry_allocate(&src.g, cfg.grid))
        a3_fail(a, "static import geometry allocation");
    int nx = cfg.grid.nx, ny = cfg.grid.ny, ns = cfg.grid.ns;
    for (int j = 0; j <= ny; j++)
        for (int i = 0; i <= nx; i++)
            memcpy(src.g.node + ((size_t)j * (nx + 1) + i) * (ns + 1),
                   g->node + ((size_t)(2 * j) * (2 * nx + 1) + 2 * i) * (ns + 1),
                   (ns + 1) * sizeof(V3));
    if (tube_geometry_rows(&src.g, &src.field, a->rank * ny / a->size,
                           (a->rank + 1) * ny / a->size))
        a3_fail(a, "original static import volumes");
    gather_geometry(a, src.g.cell, ny, (size_t)nx * ns * sizeof(TubeCell), NULL);
    int columns = nx * ny;
    src.first =
        a->rank * (columns / a->size) + (a->rank < columns % a->size ? a->rank : columns % a->size);
    src.count = columns / a->size + (a->rank < columns % a->size);
    src.n = src.count * ns;
    src.total = columns * ns;
    src.counts = a3_alloc(a, a->size, sizeof(int));
    src.offsets = a3_alloc(a, a->size, sizeof(int));
    for (int r = 0; r < a->size; r++) {
        src.counts[r] = (columns / a->size + (r < columns % a->size)) * ns;
        src.offsets[r] =
            (r * (columns / a->size) + (r < columns % a->size ? r : columns % a->size)) * ns;
    }
    src.rho = a3_alloc(a, src.total, sizeof(double));
    src.momentum = a3_alloc(a, src.total, sizeof(double));
    src.phi = a3_alloc(a, (size_t)(nx + 1) * (ny + 1), sizeof(double));
    /* Validate the full checkpoint signature; no static transport is done. */
    src.c.exb_pdm4 = a->c.exb_pdm4;
    solver3_restart(&src, a->c.restart, time, step);
    solver3_gather(&src, src.rho);
    solver3_gather(&src, src.momentum);
    double source_mass = solver3_mass(&src, src.rho);
    Config3 transfer = a->c;
    transfer.exb_pdm4 = 0;
    A3 coarse;
    a3_mesh_mask(&coarse, g, transfer, a->rank, a->size, 0, NULL, a, NULL);
    coarse.material_min = INFINITY;
    coarse.material_max = -INFINITY;
    coarse.boundary_total = src.boundary_total;
    double *N = a3_alloc(&coarse, (size_t)coarse.slots * coarse.ns, sizeof(double));
    double *M = a3_alloc(&coarse, (size_t)coarse.slots * coarse.ns, sizeof(double));
    double adjustment = 0;
    for (int slot = 0; slot < coarse.owned_count; slot++)
        for (int k = 0; k < coarse.ns; k++) {
            int id = coarse.global[slot] * coarse.ns + k, t = slot * coarse.ns + k;
            double scale = src.g.cell[id].volume / coarse.metric[id].volume;
            N[t] = src.rho[id] * scale;
            M[t] = src.momentum[id] * scale;
            adjustment = fmax(adjustment, fabs(scale - 1));
        }
    a3_remap(&coarse, a, N, M, n, m);
    a->import_mass_error = fabs(a3_mass(a, n) - source_mass) / source_mass;
    MPI_Allreduce(MPI_IN_PLACE, &adjustment, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    a->geometry_density_adjustment = adjustment;
    a->max_residual = src.max_potential_residual;
    *mass0 = src.mass0;
    if (a->import_mass_error > 1e-11)
        a3_fail(a, "static-to-AMR particle transfer failed");
    if (!a->rank) {
        FILE *f = output(a, "import_info.json", "w");
        fprintf(
            f,
            "{\"time_s\":%.17g,\"step\":%d,\"source_mass\":%.17g,\"original_mass\":%.17g,\"boundary_total\":%.17g,\"relative_particle_transfer_error\":%.17g,\"max_geometry_density_adjustment\":%.17g,\"potential_recomputed\":true}\n",
            *time, *step, source_mass, *mass0, a->boundary_total, a->import_mass_error, adjustment);
        fclose(f);
    }
    free(N);
    free(M);
    a3_free(&coarse);
    solver3_free(&src);
}
/**
 * Replace scratch arrays after a mesh rebuild so their sizes match the new owned-plus-halo
 * backing layout.
 */
static void workspace(A3 *a, double ***arrays, int count) {
    for (int i = 0; i < count; i++) {
        free(*arrays[i]);
        *arrays[i] = a3_alloc(a, (size_t)a->slots * a->ns, sizeof(double));
    }
}
/**
 * Drive initialization/restart, regridding, half chemistry, two-stage O+/molecular transport,
 * acceptance checks, diagnostics and checkpointing. All ranks follow the same clocks; no level
 * subcycling is used.
 */
int amr3_run(Config3 c, int rank, int size) {
    A3 status = {.rank = rank, .size = size};
    if (c.amr_levels == 3) {
        c.grid.nx *= 2;
        c.grid.ny *= 2;
    }
    if ((c.background_east != 0 || c.driver) &&
        ((c.periodic_x && c.amr_case != 2) || !c.regional_x || c.amr_case == 1 || c.amr_case == 3 ||
         !c.open_apex))
        a3_fail(
            &status,
            "external potential requires regional nonperiodic longitude and open-apex particles");
    if ((c.amr_case == 1 || c.amr_case == 2) && c.collision_profile[0])
        a3_fail(&status, "synthetic verification geometry has no physical collision profile");
    if ((size_t)c.grid.nx * c.grid.ny * c.grid.ns > 2000000)
        a3_fail(&status, "canonical fine geometry exceeds the current AMR size limit");
    if (c.amr_case == 1 || c.amr_case == 2) {
        c.parallel = 0;
        c.periodic_x = c.amr_case == 2 && (c.background_east != 0 || c.driver) ? 0 : 1;
    }
    if (c.amr_case == 2 || c.amr_case == 3)
        c.end = 0;
    TubeGrid g;
    geometry(&status, &g, c);
    if (c.geometry_only) {
        tube_free(&g);
        return 0;
    }
    if (c.spatial_neutral) {
        double *bath = NULL;
        if (neutral_grid_read(c.neutral_profile, &g, c.chem_Ti, c.chem_Te, &bath) != 1)
            a3_fail(&status, "spatial neutral grid/header/positions/temperatures/checksum");
        tube_conductivity_values(&g, c.e_peak, bath);
        free(bath);
        if (!rank) {
            fprintf(stdout, "Applied3D neutral bath before AMR metrics and current responses\n");
            fflush(stdout);
        }
    }
    A3 a;
    double *n = NULL, *m = NULL, mass0 = 0, time = 0, next_regrid = 0, travel = 0;
    int step = 0, retries = 0, *quiet = NULL;
    int native = a3_checkpoint_read(&a, &g, c, rank, size, &n, &m, &time, &step, &mass0, &quiet,
                                    &next_regrid, &travel, &retries);
    if (c.amr_case && c.restart[0] && !native)
        a3_fail(&status, "physical static-state import requires amr-case esf");
    if (!native)
        a3_mesh(&a, &g, c, rank, size, c.amr_mode);
    if (c.amr_case == 3) {
        remap_tests(&a, &g);
        a3_free(&a);
        tube_free(&g);
        return 0;
    }
    if (!native && c.amr_case != 1)
        a3_potential_create(&a);
    size_t total = (size_t)a.slots * a.ns;
    if (!native) {
        n = a3_alloc(&a, total, sizeof(double));
        m = a3_alloc(&a, total, sizeof(double));
    }
    double *old_n = a3_alloc(&a, total, sizeof(double)),
           *old_m = a3_alloc(&a, total, sizeof(double));
    double *stage_n = a3_alloc(&a, total, sizeof(double)),
           *stage_m = a3_alloc(&a, total, sizeof(double));
    double *dn = a3_alloc(&a, total, sizeof(double)), *dm = a3_alloc(&a, total, sizeof(double)),
           *first_m = a3_alloc(&a, total, sizeof(double));
    double **scratch[] = {&old_n, &old_m, &stage_n, &stage_m, &dn, &dm, &first_m};
    double *pre_n = NULL, *pre_m = NULL, *pre_p = NULL;
    double **pre_scratch[] = {&pre_n, &pre_m, &pre_p};
    double constant_error = 0;
    if (c.amr_case == 1 && !native) {
        a3_initialize(&a, n, m, 1);
        for (int t = 0; t < a.owned_count * a.ns; t++)
            m[t] = 0;
        rhs(&a, n, m, dn, dm, time);
        for (int t = 0; t < a.owned_count * a.ns; t++)
            constant_error = fmax(constant_error, fabs(dn[t]));
        MPI_Allreduce(MPI_IN_PLACE, &constant_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (constant_error > 1e-10)
            a3_fail(&a, "constant tracer changed at a coarse-fine interface");
    }
    if (!native) {
        a3_initialize(&a, n, m, 0);
        mass0 = a3_mass(&a, n);
        if (c.restart[0])
            import_static(&a, &g, n, m, &time, &step, &mass0);
        quiet = a3_alloc(&a, a.mask_count, sizeof(int));
        next_regrid = time + c.amr_interval;
    }
    double next = 0, start = 0, compute = 0;
    int frame = 0, own = a.owned_count * a.ns;
    double start_time = time;
    int start_step = step;
    int start_retries = retries;
    next = time;
    unsigned char *mask = a3_alloc(&a, a.mask_count, 1);
    if (c.amr_structure && c.amr_seed_guard && c.amplitude != 0 && !native) {
        if (c.restart[0])
            a3_fail(
                &a,
                "a static checkpoint has no advected seed protection: use fresh initialization or explicitly disable seed guard for a late continuation");
        if (!(c.amr_seed_lambda > 0))
            a3_fail(&a, "declare the shortest protected seed wavelength");
        a.marker = a3_alloc(&a, (size_t)a.slots * a.ns, sizeof(double));
        a3_initialize_marker(&a, n);
    }
    if (c.amr_structure) {
        a.tag_time = time;
        a.tag_start = time;
    }
    double *marker = a.marker, *old_p = a3_alloc(&a, (size_t)a.slots * a.ns, sizeof(double));
    double *stage_p = a3_alloc(&a, (size_t)a.slots * a.ns, sizeof(double)),
           *dp = a3_alloc(&a, (size_t)a.slots * a.ns, sizeof(double));
    double **marker_scratch[] = {&old_p, &stage_p, &dp};
    a.marker_rhs = marker ? dp : NULL;
    a.time = time;
    a3_night_create(&a, time);
    if (c.amr_mode == 3 && !native) {
        a3_tags(&a, n, m, mask, quiet);
        if (a3_adapt(&a, &n, &m, mask, time, step, !c.restart[0])) {
            workspace(&a, scratch, 7);
            workspace(&a, marker_scratch, 3);
        }
        marker = a.marker;
        a.marker_rhs = marker ? dp : NULL;
        own = a.owned_count * a.ns;
    }
    if (a.night)
        workspace(&a, pre_scratch, 3);
    if (c.amr_case == 1 && c.amr_structure) {
        for (int p = 0; p < a.owned_count; p++)
            for (int k = 0; k < a.ns; k++) {
                int t = p * a.ns + k;
                old_n[t] = a3_active(&a, a.global[p], k) ? 1. : 0.;
                old_m[t] = 0;
            }
        rhs(&a, old_n, old_m, dn, dm, time);
        for (int t = 0; t < a.owned_count * a.ns; t++)
            constant_error = fmax(constant_error, fabs(dn[t]));
        MPI_Allreduce(MPI_IN_PLACE, &constant_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (constant_error > 1e-10)
            a3_fail(&a, "constant changed on the actual local segmented mesh");
        rhs(&a, n, m, dn, dm, time);
    }
    for (;;) {
        a.time = time;
        if (c.amr_mode == 3 && frame > 0 && time < c.end - 1e-9 &&
            (time >= next_regrid - 1e-9 ||
             travel >= (c.amr_structure ? .5 : .5 * (c.amr_buffer + 1)))) {
            a.tag_time = time;
            a3_tags(&a, n, m, mask, quiet);
            if (a3_adapt(&a, &n, &m, mask, time, step, 0)) {
                workspace(&a, scratch, 7);
                workspace(&a, marker_scratch, 3);
                if (a.night)
                    workspace(&a, pre_scratch, 3);
            }
            marker = a.marker;
            a.marker_rhs = marker ? dp : NULL;
            own = a.owned_count * a.ns;
            next_regrid = time + c.amr_interval;
            travel = 0;
        }
        if (time >= next - 1e-9 || time >= c.end - 1e-9) {
            if (time >= c.end - 1e-9 && frame > 0)
                compute = wall3() - start;
            char file[80];
            snprintf(file, sizeof(file), "interfaces_rank%04d.csv", rank);
            a.audit_file = output(&a, file, "w");
            rhs(&a, n, m, dn, dm, time);
            fclose(a.audit_file);
            a.audit_file = NULL;
            snapshot(&a, n, m, frame++, time, step, mass0, quiet, next_regrid, travel, retries);
            a3_night_snapshot(&a, n, time, frame - 1, mass0);
            if (!rank && a3_driver_active(&a)) {
                FILE *f = output(&a, "electric_driver_diagnostics.csv", frame == 1 ? "w" : "a");
                if (frame == 1)
                    fputs("time_s,east_V_m,potential_xi_gradient_V,east_reference_scale_m\n", f);
                fprintf(f, "%.17g,%.17g,%.17g,%.17g\n", time, a3_driver_field(&a),
                        a3_driver_gradient(&a), a.east_scale_m);
                fclose(f);
            }
            if (frame == 1) {
                MPI_Barrier(MPI_COMM_WORLD);
                start = wall3();
                a.t_exchange = a.t_prepare = a.t_coeff = a.t_potential = a.t_transport =
                    a.t_regrid = 0;
                a.rhs_calls = a.tag_checks = a.chemistry_calls = a.timestep_attempts =
                    a.cfl_adjustments = 0;
            }
            next += c.output;
            if (time >= c.end - 1e-9)
                break;
        }
        if (c.amr_mode == 3 && (time >= next_regrid - 1e-9 ||
                                travel >= (c.amr_structure ? .5 : .5 * (c.amr_buffer + 1))))
            continue;
        if (step >= c.max_steps)
            a3_fail(&a, "AMR step limit");
        double dt = fmin(c.dtmax, .9 * c.cfl / fmax(a.rate, 1e-30));
        dt = fmin(dt, fmin(c.end - time, next - time));
        dt = fmin(dt, electric_driver_next(c.driver, time) - time);
        if (c.amr_mode == 3)
            dt = fmin(dt, next_regrid - time);
        memcpy(old_n, n, own * sizeof(double));
        memcpy(old_m, m, own * sizeof(double));
        if (marker)
            memcpy(old_p, marker, own * sizeof(double));
        if (a.night) {
            memcpy(pre_n, n, own * sizeof(double));
            memcpy(pre_m, m, own * sizeof(double));
            if (marker)
                memcpy(pre_p, marker, own * sizeof(double));
            a3_night_begin(&a);
        }
        int ok = 0;
        double boundary_increment = 0, travel_increment = 0;
        for (int attempt = 0; attempt < 24; attempt++) {
            a.timestep_attempts++;
            memcpy(n, a.night ? pre_n : old_n, own * sizeof(double));
            memcpy(m, a.night ? pre_m : old_m, own * sizeof(double));
            if (marker) {
                memcpy(marker, a.night ? pre_p : old_p, own * sizeof(double));
                a.marker = marker;
            }
            if (a.night) {
                a3_night_restore(&a);
                a3_night_step(&a, n, m, .5 * dt);
                memcpy(old_n, n, own * sizeof(double));
                memcpy(old_m, m, own * sizeof(double));
                if (marker)
                    memcpy(old_p, marker, own * sizeof(double));
            }
            rhs(&a, n, m, dn, dm, time);
            double boundary0 = a.boundary_rate, transverse0 = a.transverse_rate;
            double next_dt = fmin(dt, .9 * c.cfl / fmax(a.rate, 1e-30));
            if (a.night && next_dt < dt * (1 - 1e-12)) {
                a.cfl_adjustments++;
                dt = next_dt;
                continue;
            }
            dt = next_dt;
            a3_night_transport_begin(&a);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
            for (int p = 0; p < a.owned_count; p++)
                for (int k = 0; k < a.ns; k++) {
                    int t = p * a.ns + k;
                    double decay, p1, p2;
                    exponential_weights(a.metric[a.global[p] * a.ns + k].nu, dt, &decay, &p1, &p2);
                    a.etd_weight[0][t] = decay;
                    a.etd_weight[1][t] = p1;
                    a.etd_weight[2][t] = p2;
                    first_m[t] = dm[t];
                    stage_n[t] = n[t] + dt * dn[t];
                    stage_m[t] = decay * m[t] + dt * p1 * dm[t];
                    if (marker)
                        stage_p[t] = marker[t] + dt * dp[t];
                }
            if (marker)
                a.marker = stage_p;
            int molecular_ok = a3_night_predict(&a, dt);
            if (a3_valid(&a, stage_n, stage_m) && molecular_ok) {
                rhs(&a, stage_n, stage_m, dn, dm, time + dt);
                double boundary1 = a.boundary_rate;
                if (dt * a.rate <= c.cfl * (1 + 1e-12)) {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
                    for (int p = 0; p < a.owned_count; p++)
                        for (int k = 0; k < a.ns; k++) {
                            int t = p * a.ns + k;
                            double decay = a.etd_weight[0][t], p1 = a.etd_weight[1][t],
                                   p2 = a.etd_weight[2][t];
                            n[t] = .5 * old_n[t] + .5 * (stage_n[t] + dt * dn[t]);
                            m[t] = decay * old_m[t] + dt * ((p1 - p2) * first_m[t] + p2 * dm[t]);
                            if (marker)
                                marker[t] = .5 * old_p[t] + .5 * (stage_p[t] + dt * dp[t]);
                        }
                    if (marker)
                        a.marker = marker;
                    molecular_ok = a3_night_correct(&a, dt);
                    if (a3_valid(&a, n, m) && molecular_ok) {
                        boundary_increment = .5 * dt * (boundary0 + boundary1);
                        travel_increment = .5 * dt * (transverse0 + a.transverse_rate);
                        ok = 1;
                        break;
                    }
                }
            }
            retries++;
            dt *= .5;
        }
        if (!ok || dt < 1e-12)
            a3_fail(&a, "positive CFL-compliant AMR step unavailable");
        a3_night_accept(&a);
        a3_night_step(&a, n, m, .5 * dt);
        a.boundary_total += boundary_increment;
        travel += travel_increment;
        time += dt;
        step++;
        if (!rank && c.progress_steps && step % c.progress_steps == 0) {
            printf("AMR progress t=%.6f step=%d retries=%d wall=%.1f\n", time, step, retries,
                   wall3() - a.wall0);
            fflush(stdout);
        }
    }
    double mass_error = fabs(a3_mass(&a, n) - mass0 - a.boundary_total + a3_night_loss(&a)) / mass0,
           error = 0, weight = 0, potential_error = 0;
    if (c.amr_case == 1) {
        for (int p = 0; p < a.owned_count; p++) {
            for (int k = 0; k < a.ns; k++) {
                if (!a3_active(&a, a.global[p], k))
                    continue;
                A3Leaf l = a3_density_leaf(&a, a.global[p], k);
                int nx = a3_nx(&a, l.level), ny = a3_ny(&a, l.level);
                double sx = sin(PI / nx) / (PI / nx), sy = sin(PI / ny) / (PI / ny),
                       sz = sin(PI / a.ns) / (PI / a.ns);
                double truth = 1 + .2 * sx * sy * sz * sin(2 * PI * ((l.i + .5) / nx - .7 * time)) *
                                       cos(2 * PI * ((l.j + .5) / ny - .3 * time)) *
                                       cos(2 * PI * ((k + .5) / a.ns - .2 * time));
                double V = a.metric[a.global[p] * a.ns + k].volume;
                error += V * fabs(n[p * a.ns + k] - truth);
                weight += V * fabs(truth);
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &weight, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        error /= weight;
    }
    if (c.amr_case == 2) {
        for (int j = 0; j <= a.fy; j++)
            for (int i = 0; i <= a.fx; i++)
                potential_error =
                    fmax(potential_error, fabs(a.phi[j * (a.fx + 1) + i] - (double)j / a.fy -
                                               a3_driver_gradient(&a) * ((double)i / a.fx - .5)));
    }
    MPI_Allreduce(MPI_IN_PLACE, &compute, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (!rank) {
        FILE *f = output(&a, "completion.json", "w");
        fprintf(
            f,
            "{\"time_s\":%.17g,\"start_time_s\":%.17g,\"steps\":%d,\"segment_steps\":%d,\"retries\":%d,\"leaf_columns\":%d,\"cells\":%d,\"potential_unknowns\":%d,\"max_potential_relative_residual\":%.17g,\"relative_mass_error\":%.17g,\"translation_relative_L1\":%.17g,\"constant_tracer_rhs_max\":%.17g,\"manufactured_potential_max_error\":%.17g,\"compute_wall_s\":%.9g,\"wall_s\":%.9g,\"state_exchange_s\":%.9g,\"state_prepare_s\":%.9g,\"tensor_assembly_s\":%.9g,\"potential_solve_scatter_s\":%.9g,\"transport_s\":%.9g,\"regrid_inclusive_s\":%.9g,\"regrids\":%d,\"refined_root_events\":%d,\"coarsened_root_events\":%d,\"max_remap_relative_error\":%.17g,\"max_potential_relative_remap_change\":%.17g}\n",
            time, start_time, step, step - start_step, retries, a.nleaf, a.nleaf * a.ns, a.ndof,
            a.max_residual, mass_error, error, constant_error, potential_error, compute,
            wall3() - a.wall0, a.t_exchange, a.t_prepare, a.t_coeff, a.t_potential, a.t_transport,
            a.t_regrid, a.regrids, a.refinements, a.coarsenings, a.remap_error,
            a.max_potential_remap_change);
        fclose(f);
        f = output(&a, "performance_counters.json", "w");
        fprintf(
            f,
            "{\"scope\":\"After initial snapshot through terminal diagnostic RHS; rank0 call counts; reported_failure_retries and segment_failure_retries exclude checkpoint history; retries_total and start_retries include history\",\"rhs_calls\":%llu,\"tag_checks\":%llu,\"chemistry_calls\":%llu,\"timestep_attempts\":%llu,\"initial_CFL_adjustments\":%llu,\"reported_failure_retries\":%d,\"segment_failure_retries\":%d,\"retries_total\":%d,\"start_retries\":%d}\n",
            a.rhs_calls, a.tag_checks, a.chemistry_calls, a.timestep_attempts, a.cfl_adjustments,
            retries - start_retries, retries - start_retries, retries, start_retries);
        fclose(f);
    }
    free(marker);
    free(old_p);
    free(stage_p);
    free(dp);
    free(n);
    free(m);
    free(old_n);
    free(old_m);
    free(stage_n);
    free(stage_m);
    free(dn);
    free(dm);
    free(first_m);
    free(mask);
    free(quiet);
    free(pre_n);
    free(pre_m);
    free(pre_p);
    a3_night_free(&a);
    a3_free(&a);
    tube_free(&g);
    return 0;
}
