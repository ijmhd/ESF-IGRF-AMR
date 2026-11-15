/*
 * Third physical level. A macro cell covers four internal roots. Only the primary backing slot
 * evolves when its macro mask is coarse; volume-weighted virtual children serve reconstruction.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Third anisotropic level. Private levels -1,0,1 map to physical L0,L1,L2.
 * Coarse states have one owner; inactive backing slots never evolve.
 */
#include "amr3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
/**
 * Return the transverse x dimension for a private level; level -1 is the extra physical coarse
 * level.
 */
int a3_nx(const A3 *a, int level) {
    return level < 0 ? a->nx / 2 : a->nx << level;
}
/**
 * Return the transverse y dimension for a private level; level -1 is the extra physical coarse
 * level.
 */
int a3_ny(const A3 *a, int level) {
    return level < 0 ? a->ny / 2 : a->ny << level;
}
/**
 * Return the number of canonical finest cells across one requested-level cell: 4, 2 or 1.
 */
int a3_scale(int level) {
    return 1 << (1 - level);
}
/**
 * Map an internal root index to its containing two-by-two macro cell.
 */
int a3_macro_index(const A3 *a, int b) {
    return (b / a->nx / 2) * (a->nx / 2) + (b % a->nx) / 2;
}
/**
 * Return the unique backing column used to evolve an unrefined macro cell.
 */
int a3_macro_primary(const A3 *a, int b) {
    return a->first[(b / a->nx / 2) * 2 * a->nx + 2 * ((b % a->nx) / 2)];
}
/**
 * Repair interior pointers after copying an A3 struct; they must point into the destination's
 * own embedded storage.
 */
void a3_rebase(A3 *a) {
    a->layout = a->layout_storage + 1;
    a->weights = a->weight_storage + 1;
    a->parallel_length = a->parallel_length_storage + 1;
    a->parallel_response = a->parallel_response_storage + 1;
}
/**
 * Wrap or clamp macro coordinates according to physical versus synthetic boundary policy.
 */
static int coord(A3 *a, int x, int n, int axis) {
    if (a->c.amr_case == 1 || (!axis && a->c.periodic_x)) {
        x %= n;
        return x < 0 ? x + n : x;
    }
    return x < 0 ? 0 : x >= n ? n - 1 : x;
}
/**
 * Apply the two-argument minmod slope limiter without changing a conserved cell average.
 */
static double mm(double a, double b) {
    return a * b > 0 ? copysign(fmin(fabs(a), fabs(b)), a) : 0;
}
/**
 * Test whether this rank holds the primary backing state for a macro cell.
 */
static int present(A3 *a, int i, int j) {
    return a->slot[a->first[(2 * j) * a->nx + 2 * i]] >= 0;
}
/**
 * Read a neighboring macro representation, falling back to the central cell if the neighbor is
 * outside the local halo.
 */
static double around(A3 *a, int f, int i, int j, int k, int ci, int cj) {
    i = coord(a, i, a->nx / 2, 0);
    j = coord(a, j, a->ny / 2, 1);
    if (!present(a, i, j)) {
        i = ci;
        j = cj;
    }
    return a->macro_state[f][(j * (a->nx / 2) + i) * a->ns + k];
}
/**
 * Allocate macro masks/history/state tables and copy their history from the old mesh; reuse
 * immutable volume/background values on rebuild.
 */
void a3_macro_init(A3 *a, A3 *reuse, int mode) {
    if (a->c.amr_levels != 3)
        return;
    size_t cells = (size_t)a->nb / 4 * a->ns;
    int nx = a->nx / 2;
    a->macro_mask = a3_alloc(a, cells, 1);
    a->macro_next = a3_alloc(a, cells, 1);
    a->macro_field = a3_alloc(a, a->nb / 4, 1);
    a->macro_quiet = a3_alloc(a, cells, sizeof(double));
    a->macro_born = a3_alloc(a, cells, sizeof(double));
    a->macro_volume = a3_alloc(a, cells, sizeof(double));
    a->macro_bg = a3_alloc(a, cells, sizeof(double));
    for (int f = 0; f < 4; f++)
        a->macro_state[f] = a3_alloc(a, cells, sizeof(double));
    memset(a->macro_mask, 1, cells);
    memset(a->macro_next, 1, cells);
    if (!reuse && a->c.amr_case) {
        for (int j = 0; j < a->ny / 2; j++)
            for (int i = 0; i < a->nx / 2; i++)
                for (int k = 0; k < a->ns; k++) {
                    double x = (i + .5) / (a->nx / 2);
                    int on = mode != 0;
                    if (mode == 1)
                        on = x > .125 && x < .875;
                    a->macro_mask[(j * (a->nx / 2) + i) * a->ns + k] = (unsigned char)on;
                    a->macro_next[(j * (a->nx / 2) + i) * a->ns + k] = (unsigned char)on;
                }
    }
    for (size_t t = 0; t < cells; t++)
        a->macro_quiet[t] = -1;
    if (reuse) {
        memcpy(a->macro_mask, reuse->macro_next, cells);
        memcpy(a->macro_next, reuse->macro_next, cells);
        memcpy(a->macro_quiet, reuse->macro_quiet, cells * sizeof(double));
        memcpy(a->macro_born, reuse->macro_born, cells * sizeof(double));
    }
    /* Geometry and the reference background do not change on a regrid.
     * Copy immutable values without repeating their four-child reductions.
     * Keep separate ownership: a3_free still frees each mesh's tables. */
    if (reuse) {
        memcpy(a->macro_volume, reuse->macro_volume, cells * sizeof(double));
        memcpy(a->macro_bg, reuse->macro_bg, cells * sizeof(double));
        return;
    }
    for (int j = 0; j < a->ny / 2; j++)
        for (int i = 0; i < nx; i++) {
            int b = j * nx + i;
            for (int k = 0; k < a->ns; k++) {
                int t = b * a->ns + k;
                for (int v = 0; v < 2; v++)
                    for (int u = 0; u < 2; u++) {
                        int r = ((2 * j + v) * a->nx + 2 * i + u) * a->ns + k;
                        double V = a->root_volume[r];
                        a->macro_volume[t] += V;
                        a->macro_bg[t] += V * a->root_bg[r];
                    }
                a->macro_bg[t] /= a->macro_volume[t];
            }
        }
}
/**
 * Reconstruct unsegmented raw metrics from fixed canonical cells before applying any density
 * mask. Repeated invocation must not aggregate already aggregated volumes.
 */
void a3_macro_metrics(A3 *a) {
    /* Reset before each mask application, including native restart. */
    for (int p = 0; p < a->nleaf; p++) {
        A3Leaf l = a->leaf[p];
        int s = a3_scale(l.level);
        for (int k = 0; k < a->ns; k++) {
            A3Metric *m = a->metric + p * a->ns + k;
            memset(m, 0, sizeof(*m));
            for (int v = 0; v < s; v++)
                for (int u = 0; u < s; u++) {
                    int i = s * l.i + u, j = s * l.j + v, id = tube_cell_id(a->g, i, j, k);
                    TubeCell *f = a->g->cell + id;
                    double V = f->volume;
                    m->volume += V;
                    m->nu += V * f->nu;
                    m->gravity += V * f->gravity_parallel;
                    m->bg += V * a->fine_bg[id];
                    for (int z = 0; z < 2; z++)
                        m->area[z] +=
                            a->g->face[2][tube_face_id(a->g, 2, i, j, k + z)].parallel_area;
                }
            m->nu /= m->volume;
            m->gravity /= m->volume;
            m->bg /= m->volume;
        }
    }
}
/**
 * Restrict evolved owner states to macro averages, then construct bounded conservative virtual
 * internal-root states for coarse macro cells.
 */
void a3_macro_prepare(A3 *a) {
    int nx = a->nx / 2, ny = a->ny / 2;
    /* First restrict true evolved states, before overwriting virtual children. */
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            if (present(a, i, j))
                for (int k = 0; k < a->ns; k++) {
                    int t = (j * nx + i) * a->ns + k;
                    double N = 0, M = 0, P = 0;
                    for (int v = 0; v < 2; v++)
                        for (int u = 0; u < 2; u++) {
                            int r = (2 * j + v) * a->nx + 2 * i + u;
                            for (int c = 0; c < (a->refined[r] ? 4 : 1); c++) {
                                int p = a->first[r] + c, slot = a->slot[p];
                                if (slot < 0)
                                    a3_fail(a, "partial macro halo");
                                double V = a->metric[p * a->ns + k].volume;
                                N += V * a->current_n[slot * a->ns + k];
                                M += V * a->current_m[slot * a->ns + k];
                                if (a->marker)
                                    P += V * a->marker[slot * a->ns + k];
                            }
                        }
                    a->macro_state[0][t] = N / a->macro_volume[t];
                    a->macro_state[1][t] = M / a->macro_volume[t];
                    a->macro_state[2][t] = N > 0 ? M / N : 0;
                    a->macro_state[3][t] = P / a->macro_volume[t];
                }
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            if (present(a, i, j))
                for (int k = 0; k < a->ns; k++) {
                    int t = (j * nx + i) * a->ns + k;
                    if (a->macro_mask[t])
                        continue;
                    double N = a->macro_state[0][t], vel = a->macro_state[2][t],
                           bg = a->macro_bg[t], num = 0, den = 0, d[4], db[4], v[4];
                    int di[] = {-1, 1, 0, 0}, dj[] = {0, 0, -1, 1};
                    double lo = N, hi = N, vl = vel, vh = vel;
                    for (int z = 0; z < 4; z++) {
                        int x = coord(a, i + di[z], nx, 0), y = coord(a, j + dj[z], ny, 1);
                        if (!present(a, x, y)) {
                            x = i;
                            y = j;
                        }
                        int tt = (y * nx + x) * a->ns + k;
                        db[z] = a->macro_bg[tt] - bg;
                        d[z] = a->macro_state[0][tt] - N;
                        num += db[z] * d[z];
                        den += db[z] * db[z];
                        v[z] = a->macro_state[2][tt];
                        lo = fmin(lo, N + d[z]);
                        hi = fmax(hi, N + d[z]);
                        vl = fmin(vl, v[z]);
                        vh = fmax(vh, v[z]);
                    }
                    int proportional = bg > 0;
                    for (int z = 0; z < 4; z++) {
                        double neighbor = around(a, 0, i + di[z], j + dj[z], k, i, j);
                        int x = coord(a, i + di[z], nx, 0), y = coord(a, j + dj[z], ny, 1);
                        if (!present(a, x, y)) {
                            x = i;
                            y = j;
                        }
                        double left = neighbor * bg,
                               right = N * a->macro_bg[(y * nx + x) * a->ns + k];
                        if (fabs(left - right) > 64 * 2.2204460492503131e-16 *
                                                     fmax(fmax(fabs(left), fabs(right)), 1e-300))
                            proportional = 0;
                    }
                    double A = proportional ? N / bg
                               : den > 0    ? num / den
                                            : 0,
                           sx = mm(-d[0] + A * db[0], d[1] - A * db[1]),
                           sy = mm(-d[2] + A * db[2], d[3] - A * db[3]);
                    /* A zero residual background basis must not acquire a cancellation
                     * residual when daughters are many orders below the parent mean. */
                    double roundoff = 64 * 2.2204460492503131e-16 * fmax(fabs(N), fabs(A * bg));
                    if (bg > 0 && fabs(N - A * bg) <= roundoff && fabs(sx) <= roundoff &&
                        fabs(sy) <= roundoff) {
                        proportional = 1;
                        A = N / bg;
                        sx = sy = 0;
                    }
                    double cx = 0, cy = 0;
                    for (int c = 0; c < 4; c++) {
                        int r = ((2 * j + (c >> 1)) * a->nx + 2 * i + (c & 1)) * a->ns + k;
                        double V = a->root_volume[r];
                        cx += V * (c & 1 ? .25 : -.25);
                        cy += V * (c >> 1 ? .25 : -.25);
                    }
                    cx /= a->macro_volume[t];
                    cy /= a->macro_volume[t];
                    double child[4], theta = 1, particles = 0, px = 0, py = 0;
                    /* The existing background-shaped representation preserves the quiet
                     * Chapman gradient. Bounds include the represented child background. */
                    for (int c = 0; c < 4; c++) {
                        int r = ((2 * j + (c >> 1)) * a->nx + 2 * i + (c & 1)) * a->ns + k;
                        child[c] = proportional ? N * (a->root_bg[r] / bg)
                                                : A * a->root_bg[r] + (N - A * bg) +
                                                      sx * ((c & 1 ? .25 : -.25) - cx) +
                                                      sy * ((c >> 1 ? .25 : -.25) - cy);
                        double delta = child[c] - N;
                        if (child[c] < a->field_min && delta < 0)
                            theta = fmin(theta, fmax(0., (a->field_min - N) / delta));
                        if (child[c] > a->field_max && delta > 0)
                            theta = fmin(theta, fmax(0., (a->field_max - N) / delta));
                        if (child[c] <= 0 && delta < 0)
                            theta = fmin(theta, .99 * N / -delta);
                    }
                    for (int c = 0; c < 4; c++) {
                        int r = ((2 * j + (c >> 1)) * a->nx + 2 * i + (c & 1)) * a->ns + k;
                        child[c] = theta == 1 ? child[c] : N + theta * (child[c] - N);
                        double W = child[c] * a->root_volume[r];
                        particles += W;
                        px += W * (c & 1 ? .25 : -.25);
                        py += W * (c >> 1 ? .25 : -.25);
                    }
                    if (particles > 0) {
                        px /= particles;
                        py /= particles;
                    }
                    double vx = mm(vel - v[0], v[1] - vel), vy = mm(vel - v[2], v[3] - vel), dv[4],
                           limit = 1;
                    for (int c = 0; c < 4; c++) {
                        dv[c] =
                            vx * ((c & 1 ? .25 : -.25) - px) + vy * ((c >> 1 ? .25 : -.25) - py);
                        if (dv[c] > 0)
                            limit = fmin(limit, (vh - vel) / dv[c]);
                        if (dv[c] < 0)
                            limit = fmin(limit, (vl - vel) / dv[c]);
                    }
                    for (int c = 0; c < 4; c++) {
                        int r = ((2 * j + (c >> 1)) * a->nx + 2 * i + (c & 1)) * a->ns + k;
                        double u = vel + limit * dv[c];
                        a->base_n[r] = child[c];
                        a->base_m[r] = child[c] * u;
                        a->base_velocity[r] = u;
                        a->base_marker[r] = N > 0 ? child[c] * a->macro_state[3][t] / N : 0;
                    }
                }
}
/**
 * Read a prepared macro density/momentum/velocity/marker sample and verify local halo
 * availability.
 */
double a3_macro_sample(A3 *a, int field, int i, int j, int k) {
    i = coord(a, i, a->nx / 2, 0);
    j = coord(a, j, a->ny / 2, 1);
    k = coord(a, k, a->ns, 2);
    if (!present(a, i, j))
        a3_fail(a, "macro stencil exceeds exchanged halo");
    return a->macro_state[field][(j * (a->nx / 2) + i) * a->ns + k];
}
/**
 * Choose coarse-level refinement from residual gradient/curvature and parent information loss;
 * enforce nesting, buffering and physical-time quiet hysteresis.
 */
void a3_macro_tags(A3 *a, const unsigned char *fine_mask) {
    int nx = a->nx / 2, ny = a->ny / 2, slab = a->c.amr_slab;
    size_t cells = (size_t)nx * ny * a->ns;
    unsigned char *want = a3_alloc(a, cells, 1), *core = a3_alloc(a, cells, 1);
    for (int j = a->row_first / 2; j < a->row_last / 2; j++)
        for (int i = 0; i < nx; i++)
            for (int k = 0; k < a->ns; k++) {
                int t = (j * nx + i) * a->ns + k;
                double N = a->macro_state[0][t], norm = fmax(1e-3 * a->density_scale, N), d[4],
                       db[4], num = 0, den = 0;
                int di[] = {-1, 1, 0, 0}, dj[] = {0, 0, -1, 1};
                for (int z = 0; z < 4; z++) {
                    int x = coord(a, i + di[z], nx, 0), y = coord(a, j + dj[z], ny, 1);
                    double q = around(a, 0, x, y, k, i, j);
                    d[z] = q - N;
                    norm = fmax(norm, q);
                    db[z] = a->macro_bg[(y * nx + x) * a->ns + k] - a->macro_bg[t];
                    num += db[z] * d[z];
                    den += db[z] * db[z];
                }
                double A = a->c.amr_residual && den > 0 ? num / den : 0;
                for (int z = 0; z < 4; z++)
                    d[z] -= A * db[z];
                double G = fmax(fabs(.5 * (d[1] - d[0])), fabs(.5 * (d[3] - d[2]))) / norm,
                       C = hypot(d[1] + d[0], d[3] + d[2]) / norm;
                int old = a->macro_mask[t],
                    need = G > (old ? a->c.amr_grad_off : a->c.amr_grad_on) ||
                           C > (old ? a->c.amr_curve_off : a->c.amr_curve_on);
                if (old) {
                    double loss = 0;
                    for (int c = 0; c < 4; c++) {
                        int r = ((2 * j + (c >> 1)) * a->nx + 2 * i + (c & 1)) * a->ns + k;
                        loss = fmax(
                            loss,
                            fabs((a->base_n[r] - N) - A * (a->root_bg[r] - a->macro_bg[t])) / norm);
                    }
                    if (loss > .005)
                        need = 1;
                }
                core[t] = (unsigned char)need;
            }
    MPI_Allreduce(MPI_IN_PLACE, core, (int)cells, MPI_UNSIGNED_CHAR, MPI_MAX, MPI_COMM_WORLD);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            for (int k = 0; k < a->ns; k++)
                if (core[(j * nx + i) * a->ns + k])
                    for (int v = -a->c.amr_buffer; v <= a->c.amr_buffer; v++)
                        for (int u = -a->c.amr_buffer; u <= a->c.amr_buffer; u++) {
                            int x = i + u, y = j + v;
                            if (a->c.periodic_x || a->c.amr_case == 1)
                                x = coord(a, x, nx, 0);
                            if (a->c.amr_case == 1)
                                y = coord(a, y, ny, 1);
                            if (x >= 0 && x < nx && y >= 0 && y < ny)
                                want[(y * nx + x) * a->ns + k] = 1;
                        }
    /* Proper nesting, including along-field slab contacts and corners. */
    for (int j = 0; j < a->ny; j++)
        for (int i = 0; i < a->nx; i++)
            for (int k = 0; k < a->ns; k++)
                if (fine_mask[(j * a->nx + i) * a->ns + k]) {
                    for (int v = -1; v <= 1; v++)
                        for (int u = -1; u <= 1; u++)
                            for (int z = -1; z <= 1; z++) {
                                int x = i / 2 + u, y = j / 2 + v, kk = k + z * slab;
                                if (a->c.amr_case == 1 || a->c.periodic_x)
                                    x = coord(a, x, nx, 0);
                                if (a->c.amr_case == 1) {
                                    y = coord(a, y, ny, 1);
                                    kk = coord(a, kk, a->ns, 2);
                                }
                                if (x >= 0 && x < nx && y >= 0 && y < ny && kk >= 0 && kk < a->ns)
                                    want[(y * nx + x) * a->ns + kk] = 1;
                            }
                }
    for (int b = 0; b < nx * ny; b++)
        for (int k0 = 0; k0 < a->ns; k0 += slab) {
            int request = 0;
            for (int k = k0; k < k0 + slab; k++)
                request |= want[b * a->ns + k];
            for (int k = k0; k < k0 + slab; k++) {
                int t = b * a->ns + k;
                if (request) {
                    a->macro_next[t] = 1;
                    a->macro_quiet[t] = -1;
                    if (!a->macro_mask[t])
                        a->macro_born[t] = a->tag_time;
                } else if (a->macro_mask[t]) {
                    if (a->macro_quiet[t] < 0)
                        a->macro_quiet[t] = a->tag_time;
                    a->macro_next[t] =
                        (unsigned char)(a->tag_time - a->macro_quiet[t] < a->c.amr_quiet_s ||
                                        a->tag_time - a->macro_born[t] < a->c.amr_quiet_s);
                } else {
                    a->macro_next[t] = 0;
                    a->macro_quiet[t] = -1;
                }
            }
        }
    free(want);
    free(core);
}
/**
 * Write physical L0/L1/L2 cell labels and actual activity counts for visualization and
 * validation.
 */
void a3_macro_output(A3 *a, int frame, double time) {
    if (a->c.amr_levels != 3 || a->rank)
        return;
    char path[1200];
    snprintf(path, sizeof(path), "%s/amr_level_%04d.bin", a->c.out, frame);
    FILE *f = fopen(path, "wb");
    if (!f)
        a3_fail(a, "three-level output open");
    unsigned long count[3] = {0};
    for (int j = 0; j < a->fy; j++)
        for (int i = 0; i < a->fx; i++)
            for (int k = 0; k < a->ns; k++) {
                int p = a->cover[j * a->fx + i];
                unsigned char lev = (unsigned char)(a3_density_leaf(a, p, k).level + 1);
                if (fwrite(&lev, 1, 1, f) != 1)
                    a3_fail(a, "three-level output write");
                count[lev]++;
            }
    fclose(f);
    snprintf(path, sizeof(path), "%s/three_level_status.json", a->c.out);
    f = fopen(path, "w");
    if (!f)
        a3_fail(a, "level status open");
    fprintf(
        f,
        "{\"time_s\":%.17g,\"physical_levels\":3,\"logical_root_nx\":%d,\"logical_root_ny\":%d,\"level_canonical_samples\":[%lu,%lu,%lu],\"active_density_cells\":%d,\"canonical_backing_cells\":%d,\"level_format\":\"uint8[fine_ny,fine_nx,ns], physical levels0,1,2\"}\n",
        time, a->nx / 2, a->ny / 2, count[0], count[1], count[2], a->active_cells,
        a->nleaf * a->ns);
    fclose(f);
}
/**
 * Flatten a transverse potential node index on the canonical nodal grid.
 */
static int node(A3 *a, int i, int j) {
    if (a->c.periodic_x && i == a->fx)
        i = 0;
    return j * (a->fx + 1) + i;
}
/**
 * Determine a macro-edge constraint in the transverse potential approximation space.
 */
static int edge(A3 *a, int i, int j) {
    return j == 0 || j == a->fy || (!a->c.periodic_x && !a->c.regional_x && (i == 0 || i == a->fx));
}
/**
 * Accumulate constrained basis weights for one potential map without losing contributions shared
 * through a parent.
 */
static void add_map(A3 *a, A3Map *m, int id, double w) {
    if (fabs(w) < 1e-15)
        return;
    for (int z = 0; z < m->n; z++)
        if (m->id[z] == id) {
            m->w[z] += w;
            return;
        }
    if (m->n == 4)
        a3_fail(a, "hanging vertex support exceeds balanced Q1 stencil");
    m->id[m->n] = id;
    m->w[m->n++] = w;
}
/**
 * Resolve nested potential-node constraints to the final independent DOFs.
 */
static void resolve(A3 *a, int id, int *origin, unsigned char *scale, unsigned char *done) {
    if (done[id] == 2)
        return;
    if (done[id] == 1)
        a3_fail(a, "cyclic hanging-node constraint");
    done[id] = 1;
    int x = id % (a->fx + 1), y = id / (a->fx + 1), b = origin[id], s = scale[id];
    if (s && !edge(a, x, y)) {
        int bx = b % (a->fx + 1), by = b / (a->fx + 1), xx = x;
        if (a->c.periodic_x && xx == 0 && bx == a->fx - s)
            xx = a->fx;
        double u = (double)(xx - bx) / s, v = (double)(y - by) / s;
        for (int c = 0; c < 4; c++) {
            double w = (c & 1 ? u : 1 - u) * (c >> 1 ? v : 1 - v);
            if (w == 0)
                continue;
            int n = node(a, bx + s * (c & 1), by + s * (c >> 1));
            resolve(a, n, origin, scale, done);
            A3Map *q = a->map + n;
            a->map[id].fixed += w * q->fixed;
            a->map[id].drive += w * q->drive;
            for (int z = 0; z < q->n; z++)
                add_map(a, a->map + id, q->id[z], w * q->w[z]);
        }
    }
    done[id] = 2;
}
/**
 * Build nested bilinear potential constraints for the extra coarse level, keeping boundary and
 * hanging-node consistency.
 */
void a3_macro_maps(A3 *a) {
    int nv = (a->fx + 1) * (a->fy + 1), nx = a->nx / 2, ny = a->ny / 2;
    int *origin = a3_alloc(a, nv, sizeof(int));
    unsigned char *scale = a3_alloc(a, nv, 1), *done = a3_alloc(a, nv, 1);
    /* Columns of the elliptic space may remain finer than local density. */
    memset(a->macro_field, 0, a->nb / 4);
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            int b = j * nx + i;
            for (int k = 0; k < a->ns; k++)
                a->macro_field[b] |= a->macro_mask[b * a->ns + k];
            for (int v = 0; v < 2; v++)
                for (int u = 0; u < 2; u++)
                    if (a->refined[(2 * j + v) * a->nx + 2 * i + u])
                        for (int yy = j - 1; yy <= j + 1; yy++)
                            for (int xx = i - 1; xx <= i + 1; xx++) {
                                int x = xx, y = yy;
                                if (a->c.periodic_x)
                                    x = coord(a, x, nx, 0);
                                if (x >= 0 && x < nx && y >= 0 && y < ny)
                                    a->macro_field[y * nx + x] = 1;
                            }
        }
    /* Mark the finest hidden vertices first, then override with coarse parents. */
    for (int j = 0; j < a->ny; j++)
        for (int i = 0; i < a->nx; i++)
            if (!a->refined[j * a->nx + i])
                for (int v = 0; v <= 2; v++)
                    for (int u = 0; u <= 2; u++)
                        if ((u & 1) || (v & 1)) {
                            int n = node(a, 2 * i + u, 2 * j + v);
                            origin[n] = node(a, 2 * i, 2 * j);
                            scale[n] = 2;
                        }
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
            if (!a->macro_field[j * nx + i])
                for (int v = 0; v <= 4; v++)
                    for (int u = 0; u <= 4; u++)
                        if (u % 4 || v % 4) {
                            int n = node(a, 4 * i + u, 4 * j + v);
                            origin[n] = node(a, 4 * i, 4 * j);
                            scale[n] = 4;
                        }
    for (int j = 0; j <= a->fy; j++)
        for (int i = 0; i < a->fx + !a->c.periodic_x; i++) {
            int n = node(a, i, j);
            if (edge(a, i, j)) {
                a->map[n].fixed = a->c.amr_case == 2 ? (double)j / a->fy : 0;
                a->map[n].drive =
                    a3_driver_active(a) ? -a->east_scale_m * ((double)i / a->fx - .5) : 0;
                done[n] = 2;
            } else if (!scale[n]) {
                add_map(a, a->map + n, a->ndof++, 1.);
                done[n] = 2;
            }
        }
    for (int j = 0; j <= a->fy; j++)
        for (int i = 0; i < a->fx + !a->c.periodic_x; i++)
            resolve(a, node(a, i, j), origin, scale, done);
    if (a->c.periodic_x)
        for (int j = 0; j <= a->fy; j++)
            a->map[j * (a->fx + 1) + a->fx] = a->map[j * (a->fx + 1)];
    free(origin);
    free(scale);
    free(done);
}
