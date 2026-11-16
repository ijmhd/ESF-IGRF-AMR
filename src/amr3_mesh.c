/*
 * Mesh, owned/ghost storage and conservative reconstruction. Backing columns are not all active
 * density DOFs. Virtual samples may lazily populate caches; this module is not generally thread-
 * safe.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only */
#include "amr3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.14159265358979323846

/**
 * Report rank and invariant violation, abort the communicator, and never continue with an
 * invalid model state.
 */
_Noreturn void a3_fail(A3 *a, const char *message) {
    fprintf(stderr, "rank %d AMR3: %s\n", a->rank, message);
    MPI_Abort(MPI_COMM_WORLD, 1);
    abort();
}
/**
 * Allocate zero-initialized storage, including a nonnull minimum allocation for empty logical
 * arrays; allocation failure is fatal.
 */
void *a3_alloc(A3 *a, size_t n, size_t width) {
    void *p = calloc(n ? n : 1, width);
    if (!p)
        a3_fail(a, "allocation failed");
    return p;
}
/**
 * Return the nonnegative periodic index i modulo n.
 */
static int wrap(int i, int n) {
    i %= n;
    return i < 0 ? i + n : i;
}
/**
 * Wrap synthetic/periodic indices and otherwise clamp physical boundary indices.
 */
static int index_axis(A3 *a, int axis, int i, int n) {
    if (a->c.amr_case == 1 || (axis == 0 && a->c.periodic_x))
        return wrap(i, n);
    return i < 0 ? 0 : i >= n ? n - 1 : i;
}
/**
 * Return the sign-consistent slope with smallest magnitude; return zero at an extremum.
 */
static double minmod(double a, double b) {
    return a * b > 0 ? copysign(fmin(fabs(a), fabs(b)), a) : 0;
}
/**
 * Flatten a density-level cell index; k is contiguous and nx depends on the private level.
 */
static int gid(A3 *a, int level, int i, int j, int k) {
    int nx = a3_nx(a, level);
    return (j * nx + i) * a->ns + k;
}
/**
 * Integrate canonical physical-volume moments, restrict them to coarser levels, and build four-
 * cell face weights for PDM4.
 */
static void moment_weights(A3 *a) {
    const int choose[4][4] = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 2, 1, 0}, {1, 3, 3, 1}};
    double *mom_storage[3], **mom = mom_storage + 1;
    mom[0] = a3_alloc(a, (size_t)a->nb * a->ns * 12, sizeof(double));
    mom[1] = a3_alloc(a, (size_t)a->fx * a->fy * a->ns * 12, sizeof(double));
    for (int j = 0; j < a->fy; j++)
        for (int i = 0; i < a->fx; i++)
            for (int k = 0; k < a->ns; k++)
                if (tube_cell_moments(a->g, i, j, k, mom[1] + 12 * gid(a, 1, i, j, k)))
                    a3_fail(a, "invalid canonical volume moments");
    for (int j = 0; j < a->ny; j++)
        for (int i = 0; i < a->nx; i++)
            for (int k = 0; k < a->ns; k++) {
                int b = (j * a->nx + i) * a->ns + k;
                double *out = mom[0] + 12 * b, V = a->root_volume[b];
                for (int c = 0; c < 4; c++) {
                    int id = gid(a, 1, 2 * i + (c & 1), 2 * j + (c >> 1), k);
                    double weight = a->g->cell[id].volume / V;
                    for (int axis = 0; axis < 3; axis++)
                        for (int p = 0; p < 4; p++) {
                            double value = 0;
                            for (int r = 0; r <= p; r++) {
                                double scale =
                                    axis == 2
                                        ? (r == p ? 1. : 0.)
                                        : choose[p][r] *
                                              pow((double)(axis == 0 ? (c & 1) : (c >> 1)), p - r) /
                                              pow(2., p);
                                value += scale * mom[1][12 * id + 4 * axis + r];
                            }
                            out[4 * axis + p] += weight * value;
                        }
                }
            }
    if (a->c.amr_levels == 3) {
        int nx = a->nx / 2, ny = a->ny / 2;
        mom[-1] = a3_alloc(a, (size_t)nx * ny * a->ns * 12, sizeof(double));
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                for (int k = 0; k < a->ns; k++) {
                    double V = 0;
                    for (int c = 0; c < 4; c++)
                        V += a->root_volume[((2 * j + (c >> 1)) * a->nx + 2 * i + (c & 1)) * a->ns +
                                            k];
                    double *out = mom[-1] + 12 * ((j * nx + i) * a->ns + k);
                    for (int c = 0; c < 4; c++) {
                        int id = ((2 * j + (c >> 1)) * a->nx + 2 * i + (c & 1)) * a->ns + k;
                        double w = a->root_volume[id] / V;
                        for (int axis = 0; axis < 3; axis++)
                            for (int p = 0; p < 4; p++)
                                for (int r = 0; r <= p; r++) {
                                    double scale =
                                        axis == 2
                                            ? (r == p ? 1. : 0.)
                                            : choose[p][r] *
                                                  pow((double)(axis == 0 ? (c & 1) : (c >> 1)),
                                                      p - r) /
                                                  pow(2., p);
                                    out[4 * axis + p] += w * scale * mom[0][12 * id + 4 * axis + r];
                                }
                    }
                }
    }
    for (int level = a->c.amr_levels == 3 ? -1 : 0; level < 2; level++) {
        TubeGrid layout = {.c = {.nx = a3_nx(a, level), .ny = a3_ny(a, level), .ns = a->ns}};
        int nx = layout.c.nx, ny = layout.c.ny, dimension[3] = {nx, ny, a->ns};
        for (int axis = 0; axis < 3; axis++) {
            int nf = tube_face_count(&layout, axis);
            a->weights[level][axis] = a3_alloc(a, 4 * (size_t)nf, sizeof(double));
            for (int f = 0; f < 4 * nf; f++)
                a->weights[level][axis][f] = NAN;
            for (int j = 0; j < ny + (axis == 1); j++)
                for (int i = 0; i < nx + (axis == 0); i++)
                    for (int k = 0; k < a->ns + (axis == 2); k++) {
                        int coordinate[3] = {i, j, k}, face = coordinate[axis];
                        int periodic = a->c.amr_case == 1 || (axis == 0 && a->c.periodic_x);
                        if (!periodic && (face < 2 || face > dimension[axis] - 2))
                            continue;
                        double translated[4][4] = {{0}};
                        for (int c = 0; c < 4; c++) {
                            int cell[3] = {i, j, k}, off = c - 2;
                            cell[axis] = index_axis(a, axis, face + off, dimension[axis]);
                            double *local = mom[level] +
                                            12 * gid(a, level, cell[0], cell[1], cell[2]) +
                                            4 * axis;
                            for (int p = 0; p < 4; p++)
                                for (int r = 0; r <= p; r++)
                                    translated[c][p] +=
                                        choose[p][r] * pow((double)off, p - r) * local[r];
                        }
                        if (tube_moment_weights((const double(*)[4])translated,
                                                a->weights[level][axis] +
                                                    4 * tube_face_id(&layout, axis, i, j, k)))
                            a3_fail(a, "singular AMR face moments");
                    }
        }
        free(mom[level]);
    }
}
/**
 * Share reference-counted immutable geometry/response/PDM tables across regrids and allocate the
 * new metric table. Three-level raw metrics are built once by a3_macro_metrics.
 */
static void geometry_cache(A3 *a, A3 *reuse) {
    size_t roots = (size_t)a->nb * a->ns;
    if (reuse) {
        a->night = reuse->night;
        a->cache_refs = reuse->cache_refs;
        (*a->cache_refs)++;
        a->root_volume = reuse->root_volume;
        a->root_bg = reuse->root_bg;
        a->fine_bg = reuse->fine_bg;
        a->offset_x = reuse->offset_x;
        a->offset_y = reuse->offset_y;
        a->response = reuse->response;
        memcpy(a->weight_storage, reuse->weight_storage, sizeof(a->weight_storage));
        memcpy(a->parallel_length_storage, reuse->parallel_length_storage,
               sizeof(a->parallel_length_storage));
        memcpy(a->parallel_response_storage, reuse->parallel_response_storage,
               sizeof(a->parallel_response_storage));
    } else {
        a->cache_refs = a3_alloc(a, 1, sizeof(int));
        *a->cache_refs = 1;
        a->root_volume = a3_alloc(a, roots, sizeof(double));
        a->root_bg = a3_alloc(a, roots, sizeof(double));
        a->fine_bg = a3_alloc(a, 4 * roots, sizeof(double));
        a->response = a3_alloc(a, 4 * roots * 12, sizeof(double));
        for (int j = 0; j < a->fy; j++)
            for (int i = 0; i < a->fx; i++)
                for (int k = 0; k < a->ns; k++) {
                    int id = gid(a, 1, i, j, k);
                    TubeCell *c = a->g->cell + id;
                    a->fine_bg[id] =
                        a->c.amr_case == 1 ? 1. : solver3_initial_density(&a->c, c->height);
                    double *r = a->response + 12 * id;
                    double xx = vdot(c->gx, c->gx), yy = vdot(c->gy, c->gy),
                           xy = vdot(c->gx, c->gy);
                    double skew = a->c.hall ? vdot(c->gx, vcross(c->b, c->gy)) : 0.;
                    double F = c->volume * c->pF, E = c->volume * c->pE;
                    double hF = c->volume * c->hF, hE = c->volume * c->hE;
                    r[0] = F * xx;
                    r[1] = F * xy + hF * skew;
                    r[2] = F * xy - hF * skew;
                    r[3] = F * yy;
                    r[4] = E * xx;
                    r[5] = E * xy + hE * skew;
                    r[6] = E * xy - hE * skew;
                    r[7] = E * yy;
                    r[8] = c->volume * c->cF[0];
                    r[9] = c->volume * c->cF[1];
                    r[10] = c->volume * c->cE[0];
                    r[11] = c->volume * c->cE[1];
                }
        a->offset_x = a3_alloc(a, roots, sizeof(double));
        a->offset_y = a3_alloc(a, roots, sizeof(double));
        for (int j = 0; j < a->ny; j++)
            for (int i = 0; i < a->nx; i++)
                for (int k = 0; k < a->ns; k++) {
                    int b = (j * a->nx + i) * a->ns + k;
                    for (int c = 0; c < 4; c++) {
                        TubeCell *f = a->g->cell + gid(a, 1, 2 * i + (c & 1), 2 * j + (c >> 1), k);
                        double bg = a->fine_bg[gid(a, 1, 2 * i + (c & 1), 2 * j + (c >> 1), k)];
                        a->root_volume[b] += f->volume;
                        a->root_bg[b] += f->volume * bg;
                        a->offset_x[b] += f->volume * ((c & 1) ? .25 : -.25);
                        a->offset_y[b] += f->volume * ((c >> 1) ? .25 : -.25);
                    }
                    a->root_bg[b] /= a->root_volume[b];
                    a->offset_x[b] /= a->root_volume[b];
                    a->offset_y[b] /= a->root_volume[b];
                }
        if (a->c.exb_pdm4)
            moment_weights(a);
        if (a->c.parallel)
            a3_parallel_geometry(a);
    }
    a->metric = a3_alloc(a, (size_t)a->nleaf * a->ns, sizeof(A3Metric));
    /* Three-level construction calls a3_macro_metrics immediately below in
     * a3_segment_metrics. That routine must also support repeated remasking
     * during restart, so it remains the single raw-metric reset/build path.
     * Building the identical global table here first only doubles the work. */
    if (a->c.amr_levels == 3)
        return;
    for (int p = 0; p < a->nleaf; p++) {
        A3Leaf l = a->leaf[p];
        int s = a3_scale(l.level), fine_i = l.i * s, fine_j = l.j * s;
        for (int k = 0; k < a->ns; k++) {
            A3Metric *m = a->metric + p * a->ns + k;
            for (int v = 0; v < s; v++)
                for (int u = 0; u < s; u++) {
                    int id = gid(a, 1, fine_i + u, fine_j + v, k);
                    TubeCell *f = a->g->cell + id;
                    double bg = a->fine_bg[id];
                    m->volume += f->volume;
                    m->nu += f->nu * f->volume;
                    m->gravity += f->gravity_parallel * f->volume;
                    m->bg += bg * f->volume;
                    for (int z = 0; z < 2; z++)
                        m->area[z] +=
                            a->g->face[2][tube_face_id(a->g, 2, fine_i + u, fine_j + v, k + z)]
                                .parallel_area;
                }
            m->nu /= m->volume;
            m->gravity /= m->volume;
            m->bg /= m->volume;
        }
    }
}
/**
 * Balance the existing density-work estimate plus canonical-field weight using whole transverse
 * root rows. Three-level cuts align to pairs; this one-dimensional policy has a finite rank
 * limit.
 */
void a3_partition_plan(A3 *a, const unsigned char *work_mask, int *partition) {
    partition[0] = 0;
    partition[a->size] = a->ny;
    /* Local density work follows active FV cells, not union-column backing.
     * When every tube is refined uniformly along s this reduces exactly to
     * the previous leaf-column cost. Canonical field work remains additive. */
    double *work_prefix = a3_alloc(a, a->ny + 1, sizeof(double));
    for (int j = 0; j < a->ny; j++) {
        double density_work = 0;
        for (int i = 0; i < a->nx; i++)
            for (int k = 0; k < a->ns; k++)
                density_work +=
                    1 +
                    3 * (work_mask ? work_mask : a->density_refined)[(j * a->nx + i) * a->ns + k];
        work_prefix[j + 1] = work_prefix[j] + density_work / a->ns + a->c.amr_work_weight * a->nx;
    }
    for (int r = 1; r < a->size; r++) {
        int lo = partition[r - 1] + 1, hi = a->ny - (a->size - r), best = lo;
        if (a->c.amr_mode == 3) {
            double target = (double)r * work_prefix[a->ny] / a->size;
            for (int j = lo + 1; j <= hi; j++)
                if (fabs(work_prefix[j] - target) < fabs(work_prefix[best] - target))
                    best = j;
            partition[r] = best;
        } else
            partition[r] = r * a->ny / a->size;
    }
    free(work_prefix);
    if (a->c.amr_levels == 3) {
        if (a->size > a->ny / 2)
            a3_fail(a, "more ranks than physical L0 rows");
        for (int r = 1; r < a->size; r++) {
            int q = 2 * (partition[r] / 2), lo = partition[r - 1] + 2,
                hi = a->ny - 2 * (a->size - r);
            partition[r] = q < lo ? lo : q > hi ? hi : q;
        }
    }
}
/**
 * Construct a mesh from the requested fixed/uniform/dynamic mode without an explicit supplied
 * mask.
 */
void a3_mesh(A3 *a, TubeGrid *g, Config3 c, int rank, int size, int mode) {
    a3_mesh_mask(a, g, c, rank, size, mode, NULL, NULL, NULL);
}
/**
 * Construct a mesh with explicit density masks and optional old-mesh geometry reuse / ownership
 * cuts.
 */
void a3_mesh_mask(A3 *a, TubeGrid *g, Config3 c, int rank, int size, int mode,
                  const unsigned char *mask, A3 *reuse, const int *cuts) {
    a3_mesh_mask_work(a, g, c, rank, size, mode, mask, reuse, cuts, NULL);
}
/**
 * Construct backing columns, rank ownership, halos, reconstruction workspaces, geometry metrics
 * and face topology from explicit masks. Collective layout decisions must agree on all ranks.
 */
void a3_mesh_mask_work(A3 *a, TubeGrid *g, Config3 c, int rank, int size, int mode,
                       const unsigned char *mask, A3 *reuse, const int *cuts,
                       const unsigned char *work_mask) {
    memset(a, 0, sizeof(*a));
    a->layout = a->layout_storage + 1;
    a->weights = a->weight_storage + 1;
    a->parallel_length = a->parallel_length_storage + 1;
    a->parallel_response = a->parallel_response_storage + 1;
    a->c = c;
    a->g = g;
    a->rank = rank;
    a->size = size;
    a->nx = c.grid.nx;
    a->ny = c.grid.ny;
    a->fx = 2 * a->nx;
    a->fy = 2 * a->ny;
    for (int l = c.amr_levels == 3 ? -1 : 0; l < 2; l++)
        a->layout[l].c = (GridConfig){.nx = a3_nx(a, l), .ny = a3_ny(a, l), .ns = c.grid.ns};
    a->ns = c.grid.ns;
    a->nb = a->nx * a->ny;
    a->wall0 = wall3();
    a->mask_count = a->nb * (c.amr_slab ? a->ns : 1);
    a->density_refined = a3_alloc(a, (size_t)a->nb * a->ns, 1);
    a->tag_reason = a3_alloc(a, (size_t)a->nb * a->ns, sizeof(unsigned short));
    a->field_requested = a3_alloc(a, a->nb, 1);
    a->seed_mask = a3_alloc(a, (size_t)a->nb * a->ns, 1);
    a->quiet_since = a3_alloc(a, a->mask_count, sizeof(double));
    a->born = a3_alloc(a, a->mask_count, sizeof(double));
    a->tag_values = a3_alloc(a, (size_t)a->nb * a->ns * 6, sizeof(double));
    for (int b = 0; b < a->mask_count; b++)
        a->quiet_since[b] = -1;
    a->seed_lambda = c.amr_seed_lambda;
    a->branch_lambda = c.amr_branch_lambda;
    if (reuse) {
        a->tag_time = reuse->tag_time;
        a->tag_start = reuse->tag_start;
        a->seed_lambda = reuse->seed_lambda;
        a->branch_lambda = reuse->branch_lambda;
        memcpy(a->quiet_since, reuse->quiet_since, a->mask_count * sizeof(double));
        memcpy(a->born, reuse->born, a->mask_count * sizeof(double));
        memcpy(a->tag_reason, reuse->tag_reason, (size_t)a->nb * a->ns * sizeof(unsigned short));
        memcpy(a->field_requested, reuse->field_requested, a->nb);
        memcpy(a->tag_values, reuse->tag_values, (size_t)a->nb * a->ns * 6 * sizeof(double));
        memcpy(a->seed_mask, reuse->seed_mask, (size_t)a->nb * a->ns);
    }

    if (reuse) {
        a->east_scale_m = reuse->east_scale_m;
        a->time = reuse->time;
    } else if (c.amr_case == 2)
        a->east_scale_m = 1;
    else if (c.seed_east > 0 || c.background_east != 0 || c.driver) {
        Igrf field;
        V3 ref, xl, xr, yl, yr;
        double dl = 1e-4, dh = 10;
        if (solver3_load_field(&c, &field) ||
            tube_apex_label(&field, c.grid.lon + c.seed_offset_deg, c.seed_height, &ref) ||
            tube_apex_label(&field, c.grid.lon + c.seed_offset_deg - dl, c.seed_height, &xl) ||
            tube_apex_label(&field, c.grid.lon + c.seed_offset_deg + dl, c.seed_height, &xr) ||
            tube_apex_label(&field, c.grid.lon + c.seed_offset_deg, c.seed_height - dh, &yl) ||
            tube_apex_label(&field, c.grid.lon + c.seed_offset_deg, c.seed_height + dh, &yr))
            a3_fail(a, "seed reference mapping");
        V3 ex = vmul(vsub(xr, xl), c.grid.span / (2 * dl)),
           ey = vmul(vsub(yr, yl), (c.grid.hmax - c.grid.hmin) / (2 * dh));
        V3 b = vunit(igrf_field(&field, ref)),
           gx = vmul(vcross(ey, b), 1 / vdot(ex, vcross(ey, b)));
        double lon = atan2(ref.y, ref.x), projection = vdot(gx, (V3){-sin(lon), cos(lon), 0});
        if (projection <= 0)
            a3_fail(a, "seed east reference direction");
        a->east_scale_m = 1 / projection;
    }
    if (size > a->ny)
        a3_fail(a, "more ranks than root rows");
    a->refined = a3_alloc(a, a->nb, 1);
    a->first = a3_alloc(a, a->nb, sizeof(int));
    a->cover = a3_alloc(a, a->fx * a->fy, sizeof(int));
    a->row = a3_alloc(a, a->ny + 1, sizeof(int));
    a->leaf = a3_alloc(a, 4 * a->nb, sizeof(A3Leaf));
    for (int j = 0; j < a->ny; j++) {
        a->row[j] = a->nleaf;
        for (int i = 0; i < a->nx; i++) {
            int b = j * a->nx + i;
            double x = (i + .5) / a->nx, y = (j + .5) / a->ny;
            int split = 0;
            for (int k = 0; k < a->ns; k++) {
                int local =
                    mask ? mask[c.amr_slab ? b * a->ns + k : b]
                         : (mode == 2 || (mode == 1 && x > .25 && x < .75 && y > .10 && y < .65));
                a->density_refined[b * a->ns + k] = (unsigned char)local;
                split |= local;
                a->active_cells += local ? 4 : 1;
            }
            if (c.amr_structure)
                split |= a->field_requested[b];
            a->refined[b] = (unsigned char)split;
            a->split += split;
            a->first[b] = a->nleaf;
            for (int v = 0; v < (split ? 2 : 1); v++)
                for (int u = 0; u < (split ? 2 : 1); u++) {
                    int p = a->nleaf++;
                    a->leaf[p] = (A3Leaf){split, split ? 2 * i + u : i, split ? 2 * j + v : j, b};
                    for (int t = 0; t < (split ? 1 : 2); t++)
                        for (int s = 0; s < (split ? 1 : 2); s++)
                            a->cover[(2 * j + v + t) * a->fx + 2 * i + u + s] = p;
                }
        }
    }
    a->row[a->ny] = a->nleaf;
    a->partition = a3_alloc(a, size + 1, sizeof(int));
    a->partition[size] = a->ny;
    if (cuts)
        memcpy(a->partition, cuts, (size + 1) * sizeof(int));
    else
        a3_partition_plan(a, work_mask, a->partition);
    a->row_first = a->partition[rank];
    a->row_last = a->partition[rank + 1];
    a->owned_first = a->row[a->row_first];
    a->owned_count = a->row[a->row_last] - a->owned_first;
    a->slot = a3_alloc(a, a->nleaf, sizeof(int));
    a->global = a3_alloc(a, a->nleaf, sizeof(int));
    for (int p = 0; p < a->nleaf; p++)
        a->slot[p] = -1;
    for (int p = a->owned_first; p < a->owned_first + a->owned_count; p++) {
        a->slot[p] = a->slots;
        a->global[a->slots++] = p;
    }
    a->halo = a3_alloc(a, size, sizeof(A3Halo));
    a->requests = a3_alloc(a, 2 * size, sizeof(MPI_Request));
    a->statuses = a3_alloc(a, 2 * size, sizeof(MPI_Status));
    for (int peer = 0; peer < size; peer++)
        if (peer != rank) {
            int halo = c.amr_levels == 3 ? 6 : 3;
            int lo = a->partition[peer], hi = a->partition[peer + 1];
            int first = lo > a->row_first - halo ? lo : a->row_first - halo,
                last = hi < a->row_last + halo ? hi : a->row_last + halo;
            int slo = a->row_first > lo - halo ? a->row_first : lo - halo,
                shi = a->row_last < hi + halo ? a->row_last : hi + halo;
            int recv = last > first ? a->row[last] - a->row[first] : 0,
                send = shi > slo ? a->row[shi] - a->row[slo] : 0;
            /* Synthetic periodic y needs both wrap halos, which can be disjoint.
             * Exchange whole peer slabs in that verification case only. */
            if (c.amr_case == 1) {
                first = lo;
                last = hi;
                slo = a->row_first;
                shi = a->row_last;
                recv = a->row[hi] - a->row[lo];
                send = a->owned_count;
            }
            if (!recv && !send)
                continue;
            A3Halo *h = a->halo + a->peers++;
            h->peer = peer;
            h->recv_first = recv ? a->row[first] : 0;
            h->recv_count = recv;
            h->send_first = send ? a->row[slo] : 0;
            h->send_count = send;
            h->recv = a3_alloc(a, 3 * (size_t)recv * a->ns, sizeof(double));
            h->send = a3_alloc(a, 3 * (size_t)send * a->ns, sizeof(double));
            for (int p = h->recv_first; p < h->recv_first + recv; p++) {
                a->slot[p] = a->slots;
                a->global[a->slots++] = p;
            }
        }
    size_t n = (size_t)a->nb * a->ns;
    double **array[] = {&a->base_n, &a->base_m, &a->sx_n, &a->sy_n, &a->amp_n};
    for (int t = 0; t < 5; t++)
        *array[t] = a3_alloc(a, n, sizeof(double));
    a->ready = a3_alloc(a, n, 1);
    a->base_velocity = a3_alloc(a, n, sizeof(double));
    a->velocity = a3_alloc(a, (size_t)a->slots * a->ns, sizeof(double));
    a->virtual_offset = a3_alloc(a, a->nb, sizeof(int));
    int virtual_cells = 0;
    for (int b = 0; b < a->nb; b++) {
        a->virtual_offset[b] = -1;
        if (a->slot[a->first[b]] >= 0) {
            a->virtual_offset[b] = virtual_cells;
            virtual_cells += 4 * a->ns;
        }
    }
    a->child_n = a3_alloc(a, virtual_cells, sizeof(double));
    a->child_m = a3_alloc(a, virtual_cells, sizeof(double));
    a->child_velocity = a3_alloc(a, virtual_cells, sizeof(double));
    a->rates = a3_alloc(a, (size_t)a->owned_count * a->ns, sizeof(double));
    for (int f = 0; f < 3; f++)
        a->etd_weight[f] = a3_alloc(a, (size_t)a->owned_count * a->ns, sizeof(double));
    a->transverse_rates = a3_alloc(a, (size_t)a->owned_count * a->ns, sizeof(double));
    if (c.parallel) {
        size_t faces = (size_t)a->owned_count * (a->ns + 1);
        a->parallel_n = a3_alloc(a, faces, sizeof(double));
        a->parallel_m = a3_alloc(a, faces, sizeof(double));
        a->parallel_rate = a3_alloc(a, faces, sizeof(double));
    }
    geometry_cache(a, reuse);
    if (c.amr_case != 1)
        for (int d = 0; d < 3; d++)
            a->potential_gradient[d] =
                a3_alloc(a, (size_t)2 * (a->fx + 1) * (a->fy + 1), sizeof(double));
    a->density_scale = 1;
    if (c.amr_case != 1) {
        a->density_scale = 0;
        for (int t = 0; t < a->fx * a->fy * a->ns; t++)
            a->density_scale = fmax(a->density_scale, a->fine_bg[t]);
    }
    a->base_marker = a3_alloc(a, (size_t)a->nb * a->ns, sizeof(double));
    a3_macro_init(a, reuse, mode);
    a3_segment_metrics(a);
    a3_face_topology(a);
    if (a->night)
        a3_night_refresh(a);
}
/**
 * Reset and aggregate physical metrics for the actual segmented density owners; inactive backing
 * slots receive zero volume.
 */
void a3_segment_metrics(A3 *a) {
    if (a->c.amr_levels == 3)
        a3_macro_metrics(a);
    /* In a segmented coarse cell only its primary slot is an evolved DOF.
     * Other backing slots have zero volume and never supply independent states. */
    if (a->c.amr_slab || a->c.amr_structure)
        for (int b = 0; b < a->nb; b++)
            if (a->refined[b])
                for (int k = 0; k < a->ns; k++)
                    if (!a->density_refined[b * a->ns + k]) {
                        A3Metric sum = {0};
                        int first = a->first[b];
                        for (int d = 0; d < 4; d++) {
                            A3Metric q = a->metric[(first + d) * a->ns + k];
                            sum.volume += q.volume;
                            sum.nu += q.volume * q.nu;
                            sum.gravity += q.volume * q.gravity;
                            sum.bg += q.volume * q.bg;
                            sum.area[0] += q.area[0];
                            sum.area[1] += q.area[1];
                            memset(a->metric + (first + d) * a->ns + k, 0, sizeof(A3Metric));
                        }
                        sum.nu /= sum.volume;
                        sum.gravity /= sum.volume;
                        sum.bg /= sum.volume;
                        a->metric[first * a->ns + k] = sum;
                    }

    a->active_cells = 0;
    for (int b = 0; b < a->nb * a->ns; b++)
        a->active_cells += a->density_refined[b] ? 4 : 1;
    if (a->c.amr_levels == 3) {
        for (int j = 0; j < a->ny / 2; j++)
            for (int i = 0; i < a->nx / 2; i++)
                for (int k = 0; k < a->ns; k++)
                    if (!a->macro_mask[(j * (a->nx / 2) + i) * a->ns + k]) {
                        int primary = a->first[2 * j * a->nx + 2 * i];
                        A3Metric sum = {0};
                        for (int v = 0; v < 2; v++)
                            for (int u = 0; u < 2; u++) {
                                int b = (2 * j + v) * a->nx + 2 * i + u;
                                if (a->density_refined[b * a->ns + k])
                                    a3_fail(a, "L2 without parent L1");
                                int p = a->first[b];
                                A3Metric q = a->metric[p * a->ns + k];
                                sum.volume += q.volume;
                                sum.nu += q.volume * q.nu;
                                sum.gravity += q.volume * q.gravity;
                                sum.bg += q.volume * q.bg;
                                sum.area[0] += q.area[0];
                                sum.area[1] += q.area[1];
                                memset(a->metric + p * a->ns + k, 0, sizeof(A3Metric));
                            }
                        sum.nu /= sum.volume;
                        sum.gravity /= sum.volume;
                        sum.bg /= sum.volume;
                        a->metric[primary * a->ns + k] = sum;
                        a->active_cells -= 3;
                    }
    }
}

/**
 * Exchange owned O+ density/momentum and optional marker halos using nonblocking peer transfers,
 * then unpack only after completion.
 */
static void exchange(A3 *a, double *n, double *m) {
    double begin = wall3();
    int q = 0, fields = a->marker ? 3 : 2;
    for (int h = 0; h < a->peers; h++)
        if (a->halo[h].recv_count)
            MPI_Irecv(a->halo[h].recv, fields * a->halo[h].recv_count * a->ns, MPI_DOUBLE,
                      a->halo[h].peer, 131, MPI_COMM_WORLD, a->requests + q++);
    for (int h = 0; h < a->peers; h++) {
        A3Halo *p = a->halo + h;
        int count = p->send_count * a->ns;
        if (count) {
            int start = a->slot[p->send_first] * a->ns;
            memcpy(p->send, n + start, count * sizeof(double));
            memcpy(p->send + count, m + start, count * sizeof(double));
            if (a->marker)
                memcpy(p->send + 2 * count, a->marker + start, count * sizeof(double));
            MPI_Isend(p->send, fields * count, MPI_DOUBLE, p->peer, 131, MPI_COMM_WORLD,
                      a->requests + q++);
        }
    }
    MPI_Waitall(q, a->requests, a->statuses);
    for (int h = 0; h < a->peers; h++) {
        A3Halo *p = a->halo + h;
        int count = p->recv_count * a->ns;
        if (count) {
            int start = a->slot[p->recv_first] * a->ns;
            memcpy(n + start, p->recv, count * sizeof(double));
            memcpy(m + start, p->recv + count, count * sizeof(double));
            if (a->marker)
                memcpy(a->marker + start, p->recv + 2 * count, count * sizeof(double));
        }
    }
    a->t_exchange += wall3() - begin;
}
/**
 * Bind the current stage, exchange state, reset lazy caches, restrict actual conserved states,
 * construct macro representations and refresh molecular/electron responses. Must finish before
 * potential or flux sampling.
 */
void a3_prepare(A3 *a, double *n, double *m) {
    exchange(a, n, m);
    double begin = wall3();
    a->current_n = n;
    a->current_m = m;
    memset(a->ready, 0, (size_t)a->nb * a->ns);
    double low = a->material_min, high = a->material_max;
    for (int p = 0; p < a->owned_count * a->ns; p++) {
        low = fmin(low, n[p]);
        high = fmax(high, n[p]);
    }
    double bounds[2] = {-low, high};
    MPI_Allreduce(MPI_IN_PLACE, bounds, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    a->field_min = -bounds[0];
    a->field_max = bounds[1];
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int t = 0; t < a->slots * a->ns; t++) {
        a->velocity[t] = n[t] > 0 ? m[t] / n[t] : 0.;
    }
    for (int b = 0; b < a->nb; b++) {
        int first = a->first[b];
        if (a->slot[first] < 0)
            continue;
        for (int k = 0; k < a->ns; k++) {
            double N = 0, M = 0, P = 0;
            for (int c = 0; c < (a->refined[b] ? 4 : 1); c++) {
                int p = first + c, slot = a->slot[p];
                if (slot < 0)
                    a3_fail(a, "incomplete parent halo");
                double V = a->metric[p * a->ns + k].volume;
                N += V * n[slot * a->ns + k];
                M += V * m[slot * a->ns + k];
                if (a->marker)
                    P += V * a->marker[slot * a->ns + k];
            }
            int t = b * a->ns + k;
            if (a->refined[b]) {
                a->base_n[t] = N / a->root_volume[t];
                a->base_m[t] = M / a->root_volume[t];
            } else {
                int slot = a->slot[first] * a->ns + k;
                a->base_n[t] = n[slot];
                a->base_m[t] = m[slot];
            }
            a->base_marker[t] = a->marker ? P / a->root_volume[t] : 0;
            a->base_velocity[t] = a->base_n[t] > 0 ? a->base_m[t] / a->base_n[t] : 0.;
        }
    }
    if (a->c.amr_levels == 3)
        a3_macro_prepare(a);
    a->t_prepare += wall3() - begin;
    a3_night_refresh(a);
    a3_night_electrons(a);
}
/**
 * Read a parent density/momentum sample after applying boundary indices; reject a stencil
 * outside the exchanged root-row halo.
 */
static double parent(A3 *a, int field, int i, int j, int k) {
    i = index_axis(a, 0, i, a->nx);
    j = index_axis(a, 1, j, a->ny);
    k = index_axis(a, 2, k, a->ns);
    int b = j * a->nx + i;
    if (a->slot[a->first[b]] < 0)
        a3_fail(a, "state stencil exceeds the exchanged root-row halo");
    return (field ? a->base_m : a->base_n)[b * a->ns + k];
}
/**
 * Evaluate a conservative background-shaped density child from cached parent amplitude and
 * residual slopes.
 */
static double child_proposal(A3 *a, int b, int child, int k) {
    int i = b % a->nx, j = b / a->nx, t = b * a->ns + k;
    int id = gid(a, 1, 2 * i + (child & 1), 2 * j + (child >> 1), k);
    double bg = a->fine_bg[id];
    double A = a->amp_n[t];
    double sx = a->sx_n[t], sy = a->sy_n[t];
    return A * bg + ((a->ready[t] & 4) ? 0 : (a->base_n[t] - A * a->root_bg[t])) +
           sx * ((child & 1 ? .25 : -.25) - a->offset_x[t]) +
           sy * ((child >> 1 ? .25 : -.25) - a->offset_y[t]);
}
/**
 * Build bounded background-aware density reconstruction and lazy child-density caches for one
 * parent and along-field index.
 */
static void slopes(A3 *a, int b, int k) {
    int i = b % a->nx, j = b / a->nx, t = b * a->ns + k;
    int neighbor[4] = {
        j * a->nx + index_axis(a, 0, i - 1, a->nx), j * a->nx + index_axis(a, 0, i + 1, a->nx),
        index_axis(a, 1, j - 1, a->ny) * a->nx + i, index_axis(a, 1, j + 1, a->ny) * a->nx + i};
    const int di[4] = {-1, 1, 0, 0}, dj[4] = {0, 0, -1, 1};
    double bg = a->root_bg[t], delta_bg[4], difference[4], norm = 0, A = 0;
    double parent_n = parent(a, 0, i, j, k);
    for (int d = 0; d < 4; d++) {
        delta_bg[d] = a->root_bg[neighbor[d] * a->ns + k] - bg;
        difference[d] = parent(a, 0, i + di[d], j + dj[d], k) - parent_n;
        norm += delta_bg[d] * delta_bg[d];
        A += difference[d] * delta_bg[d];
    }
    A = norm > 0 ? A / norm : 0;
    if (a->c.amr_levels == 3 && bg > 0) {
        int proportional = 1;
        for (int d = 0; d < 4; d++) {
            double left = parent(a, 0, i + di[d], j + dj[d], k) * bg,
                   right = parent_n * a->root_bg[neighbor[d] * a->ns + k];
            if (fabs(left - right) >
                64 * 2.2204460492503131e-16 * fmax(fmax(fabs(left), fabs(right)), 1e-300))
                proportional = 0;
        }
        if (proportional) {
            A = parent_n / bg;
            a->ready[t] |= 4;
        }
    }
    a->amp_n[t] = A;
    a->sx_n[t] = minmod(-difference[0] + A * delta_bg[0], difference[1] - A * delta_bg[1]);
    a->sy_n[t] = minmod(-difference[2] + A * delta_bg[2], difference[3] - A * delta_bg[3]);
    if (a->ready[t] & 4)
        a->sx_n[t] = a->sy_n[t] = 0;
    double theta = 1, proposed[4];
    for (int c = 0; c < 4; c++) {
        proposed[c] = child_proposal(a, b, c, k);
        double delta = proposed[c] - parent_n;
        if (proposed[c] < a->field_min && delta < 0)
            theta = fmin(theta, fmax(0., (a->field_min - parent_n) / delta));
        if (proposed[c] > a->field_max && delta > 0)
            theta = fmin(theta, fmax(0., (a->field_max - parent_n) / delta));
        if (proposed[c] <= 0 && delta < 0)
            theta = fmin(theta, .99 * parent_n / (-delta));
    }
    int off = a->virtual_offset[b] + 4 * k;
    if (a->virtual_offset[b] < 0)
        a3_fail(a, "missing virtual-cell storage");
    for (int c = 0; c < 4; c++) {
        double N = theta == 1 ? proposed[c] : parent_n + theta * (proposed[c] - parent_n);
        a->child_n[off + c] = N;
    }
    a->ready[t] |= 1;
}
/**
 * Prolong parallel velocity with particle-weighted centering so child momentum sums to the
 * parent momentum.
 */
static void velocity_slopes(A3 *a, int b, int k) {
    int i = b % a->nx, j = b / a->nx, t = b * a->ns + k, off = a->virtual_offset[b] + 4 * k;
    if (!(a->ready[t] & 1))
        slopes(a, b, k);
    int neighbor[4] = {
        j * a->nx + index_axis(a, 0, i - 1, a->nx), j * a->nx + index_axis(a, 0, i + 1, a->nx),
        index_axis(a, 1, j - 1, a->ny) * a->nx + i, index_axis(a, 1, j + 1, a->ny) * a->nx + i};
    double particles = 0, center_x = 0, center_y = 0;
    for (int c = 0; c < 4; c++) {
        double weight = a->g->cell[gid(a, 1, 2 * i + (c & 1), 2 * j + (c >> 1), k)].volume *
                        a->child_n[off + c];
        particles += weight;
        center_x += weight * (c & 1 ? .25 : -.25);
        center_y += weight * (c >> 1 ? .25 : -.25);
    }
    double velocity = a->base_velocity[t], umin = velocity, umax = velocity, around[4];
    for (int d = 0; d < 4; d++) {
        around[d] = a->base_velocity[neighbor[d] * a->ns + k];
        umin = fmin(umin, around[d]);
        umax = fmax(umax, around[d]);
    }
    double sx = minmod(velocity - around[0], around[1] - velocity);
    double sy = minmod(velocity - around[2], around[3] - velocity);
    if (particles > 0) {
        center_x /= particles;
        center_y /= particles;
    }
    double delta[4], limit = 1;
    for (int c = 0; c < 4; c++) {
        delta[c] = sx * ((c & 1 ? .25 : -.25) - center_x) + sy * ((c >> 1 ? .25 : -.25) - center_y);
        if (delta[c] > 0)
            limit = fmin(limit, fmax(0., (umax - velocity) / delta[c]));
        if (delta[c] < 0)
            limit = fmin(limit, fmax(0., (umin - velocity) / delta[c]));
    }
    for (int c = 0; c < 4; c++) {
        a->child_velocity[off + c] = velocity + limit * delta[c];
        a->child_m[off + c] = a->child_n[off + c] * a->child_velocity[off + c];
    }
    a->ready[t] |= 2;
}

/**
 * Return a prospective child density/momentum for transfer/coarsening assessment while
 * preserving the existing reconstruction semantics.
 */
double a3_proposed_child(A3 *a, int b, int child, int k, int field) {
    int t = b * a->ns + k;
    if (!(a->ready[t] & 1))
        slopes(a, b, k);
    if (field == 3)
        return a->child_n[a->virtual_offset[b] + 4 * k + child] *
               (a->base_n[t] > 0 ? a->base_marker[t] / a->base_n[t] : 0);
    if (field && !(a->ready[t] & 2))
        velocity_slopes(a, b, k);
    return (field == 2 ? a->child_velocity
            : field    ? a->child_m
                       : a->child_n)[a->virtual_offset[b] + 4 * k + child];
}
/**
 * Sample actual or virtual states at a requested level. It may update lazy density/velocity
 * child caches; concurrent calls require explicit cache ownership.
 */
double a3_sample(A3 *a, int field, int level, int i, int j, int k) {
    if (level < 0)
        return a3_macro_sample(a, field, i, j, k);
    i = index_axis(a, 0, i, a3_nx(a, level));
    j = index_axis(a, 1, j, a3_ny(a, level));
    k = index_axis(a, 2, k, a->ns);
    if (!level) {
        if (field == 3)
            return a->base_marker[(j * a->nx + i) * a->ns + k];
        if (field == 2)
            return a->base_velocity[(j * a->nx + i) * a->ns + k];
        return parent(a, field, i, j, k);
    }
    int b = (j / 2) * a->nx + i / 2, child = (i & 1) + 2 * (j & 1);
    if (a->density_refined[b * a->ns + k]) {
        int p = a->first[b] + child, slot = a->slot[p];
        if (slot < 0)
            a3_fail(a, "missing fine halo");
        return (field == 3   ? a->marker
                : field == 2 ? a->velocity
                : field      ? a->current_m
                             : a->current_n)[slot * a->ns + k];
    }
    int t = b * a->ns + k;
    if (field == 3)
        return a3_proposed_child(a, b, child, k, field);
    if (!(a->ready[t] & 1))
        slopes(a, b, k);
    if (field && !(a->ready[t] & 2))
        velocity_slopes(a, b, k);
    return (field == 2 ? a->child_velocity
            : field    ? a->child_m
                       : a->child_n)[a->virtual_offset[b] + 4 * k + child];
}

/**
 * Collect a one-dimensional face stencil at a common requested level, including virtual samples
 * across coarse/fine contacts.
 */
static void sample_pair(A3 *a, int level, int i, int j, int k, double value[2]) {
    if (level < 0) {
        value[0] = a3_sample(a, 0, level, i, j, k);
        value[1] = a3_sample(a, a->c.parallel_relaxation ? 2 : 1, level, i, j, k);
        return;
    }
    i = index_axis(a, 0, i, a3_nx(a, level));
    j = index_axis(a, 1, j, a3_ny(a, level));
    k = index_axis(a, 2, k, a->ns);
    if (!level) {
        int b = j * a->nx + i, t = b * a->ns + k;
        if (a->slot[a->first[b]] < 0)
            a3_fail(a, "missing coarse face halo");
        value[0] = a->base_n[t];
        value[1] = (a->c.parallel_relaxation ? a->base_velocity : a->base_m)[t];
        return;
    }
    int b = (j / 2) * a->nx + i / 2, child = (i & 1) + 2 * (j & 1);
    if (a->density_refined[b * a->ns + k]) {
        int slot = a->slot[a->first[b] + child];
        if (slot < 0)
            a3_fail(a, "missing fine face halo");
        value[0] = a->current_n[slot * a->ns + k];
        value[1] = (a->c.parallel_relaxation ? a->velocity : a->current_m)[slot * a->ns + k];
    } else {
        if (!(a->ready[b * a->ns + k] & 2))
            velocity_slopes(a, b, k);
        int t = a->virtual_offset[b] + 4 * k + child;
        value[0] = a->child_n[t];
        value[1] = (a->c.parallel_relaxation ? a->child_velocity : a->child_m)[t];
    }
}
/**
 * Return particle-density and velocity face states using physical-volume PDM4 or the configured
 * MUSCL fallback.
 */
void a3_face_pair(A3 *a, int level, int axis, int i, int j, int k, int sign, double *n,
                  double *second) {
    int fid = tube_face_id(a->layout + level, axis, i, j, k);
    double *weights = a->c.exb_pdm4 ? a->weights[level][axis] + 4 * fid : NULL;
    double result[2], value[4][2];
    if (weights && isfinite(weights[0])) {
        double high[2] = {0, 0};
        for (int c = 0; c < 4; c++) {
            int cell[3] = {i, j, k};
            cell[axis] += c - 2;
            sample_pair(a, level, cell[0], cell[1], cell[2], value[c]);
            for (int f = 0; f < 2; f++)
                high[f] += weights[c] * value[c][f];
        }
        int donor = sign > 0 ? 1 : 2, before = donor - sign, after = donor + sign;
        for (int f = 0; f < 2; f++)
            result[f] = tube_pdm_limit(high[f], value[before][f], value[donor][f], value[after][f],
                                       a->c.pdm_beta);
    } else {
        int donor[3] = {i, j, k};
        donor[axis] -= sign > 0;
        for (int c = 0; c < 3; c++) {
            int cell[3] = {donor[0], donor[1], donor[2]};
            cell[axis] += c - 1;
            sample_pair(a, level, cell[0], cell[1], cell[2], value[c]);
        }
        for (int f = 0; f < 2; f++) {
            double left = value[1][f] - value[0][f], right = value[2][f] - value[1][f];
            double slope = minmod(1.3 * left, minmod(.5 * (left + right), 1.3 * right));
            result[f] = value[1][f] + .5 * sign * slope;
        }
    }
    *n = result[0];
    *second = result[1];
}
/**
 * Evaluate the idealized Chapman/seed fixture on canonical physical geometry; used for fresh
 * initialization and verification.
 */
static double initial_fine(A3 *a, int i, int j, int k, int test) {
    TubeCell *g = a->g->cell + gid(a, 1, i, j, k);
    if (test == 1 || test == 4)
        return 1.;
    if (test == 2)
        return solver3_initial_density(&a->c, g->height);
    if (test == 3)
        return i < a->fx / 2 ? 1e-8 : 1.;
    if (a->c.amr_case == 1) {
        double x = 2 * PI * (i + .5) / a->fx, y = 2 * PI * (j + .5) / a->fy,
               z = 2 * PI * (k + .5) / a->ns;
        if (a->c.amr_fixture) {
            double X = x / (2 * PI), Y = y / (2 * PI), Z = z / (2 * PI),
                   slab = exp(-pow((Z - .5) / .06, 2));
            if (a->c.amr_fixture == 1)
                return .02;
            if (a->c.amr_fixture == 2)
                return 1 - .95 * .5 * (1 + tanh((Y - .45) / .025)) * slab;
            if (a->c.amr_fixture == 3)
                return 1 - .95 * exp(-pow((X - .5) / .08, 2) - pow((Y - .55) / .18, 2)) * slab;
            if (a->c.amr_fixture == 4) {
                double branch = .5 * (1 + tanh((Y - .52) / .035)), sep = .20 * branch;
                double stem = exp(-pow((X - .5) / .045, 2));
                double fingers =
                    fmax(exp(-pow((X - .5 - sep) / .045, 2)), exp(-pow((X - .5 + sep) / .045, 2)));
                return 1 - .95 * ((1 - branch) * stem + branch * fingers) *
                               exp(-pow((Y - .55) / .32, 8)) * slab;
            }
            if (a->c.amr_fixture == 5)
                return 1 + a->c.amplitude * sin(2 * x) * exp(-pow((Y - .5) / .1, 2)) * slab;
        }
        double sx = sin(PI / a->fx) / (PI / a->fx), sy = sin(PI / a->fy) / (PI / a->fy),
               sz = sin(PI / a->ns) / (PI / a->ns);
        return 1 + .2 * sin(x) * cos(y) * cos(z) * sx * sy * sz;
    }
    double apex = a->c.grid.hmin + (a->c.grid.hmax - a->c.grid.hmin) * (j + .5) / a->fy;
    double longitudinal = cos(2 * PI * ((i + .5) / a->fx - a->c.seed_offset_deg / a->c.grid.span));
    if (a->c.seed_east > 0) {
        double u = a->east_scale_m *
                   ((i + .5) / a->fx - .5 - a->c.seed_offset_deg / a->c.grid.span) / a->c.seed_east;
        longitudinal = (2 * u * u - 1) * exp(-u * u);
    }
    double seed =
        a->c.amplitude * longitudinal * exp(-pow((apex - a->c.seed_height) / a->c.seed_width, 2));
    return solver3_initial_density(&a->c, g->height) * (1 + seed);
}
/**
 * Initialize an advected seed-protection material quantity consistently with density; this
 * tracer influences tagging, not conductivity.
 */
void a3_initialize_marker(A3 *a, const double *n) {
    (void)n;
    for (int slot = 0; slot < a->owned_count; slot++)
        for (int k = 0; k < a->ns; k++) {
            int p = a->global[slot];
            if (!a3_active(a, p, k))
                continue;
            A3Leaf l = a3_density_leaf(a, p, k);
            int s = a3_scale(l.level);
            double total = 0;
            for (int v = 0; v < s; v++)
                for (int u = 0; u < s; u++) {
                    int i = s * l.i + u, j = s * l.j + v;
                    double envelope;
                    if (a->c.amr_case == 1)
                        envelope = exp(-pow(((j + .5) / a->fy - .5) / .1, 2) -
                                       pow(((k + .5) / a->ns - .5) / .06, 2));
                    else {
                        double apex =
                            a->c.grid.hmin + (a->c.grid.hmax - a->c.grid.hmin) * (j + .5) / a->fy;
                        envelope = exp(-pow((apex - a->c.seed_height) / a->c.seed_width, 2));
                        if (a->c.seed_east > 0) {
                            double x =
                                a->east_scale_m *
                                ((i + .5) / a->fx - .5 - a->c.seed_offset_deg / a->c.grid.span) /
                                a->c.seed_east;
                            envelope *= fmin(1., (1 + 2 * x * x) * exp(-x * x));
                        }
                    }
                    double chi = envelope >= 1e-3 ? 1. : 0.;
                    total +=
                        a->g->cell[gid(a, 1, i, j, k)].volume * initial_fine(a, i, j, k, 0) * chi;
                }
            a->marker[slot * a->ns + k] = total / a->metric[p * a->ns + k].volume;
        }
}
/**
 * Restrict the canonical idealized state to actual active volumes; leave nonowner backing states
 * inactive.
 */
void a3_initialize(A3 *a, double *n, double *m, int test) {
    a->material_min = INFINITY;
    a->material_max = -INFINITY;
    for (int j = 0; j < a->fy; j++)
        for (int i = 0; i < a->fx; i++)
            for (int k = 0; k < a->ns; k++) {
                double q = initial_fine(a, i, j, k, test);
                a->material_min = fmin(a->material_min, q);
                a->material_max = fmax(a->material_max, q);
            }
    for (int slot = 0; slot < a->owned_count; slot++) {
        int p = a->global[slot];
        for (int k = 0; k < a->ns; k++) {
            if (!a3_active(a, p, k))
                continue;
            A3Leaf l = a3_density_leaf(a, p, k);
            int s = a3_scale(l.level);
            double value = 0, momentum = 0;
            for (int v = 0; v < s; v++)
                for (int u = 0; u < s; u++) {
                    int i = s * l.i + u, j = s * l.j + v;
                    double V = a->g->cell[gid(a, 1, i, j, k)].volume;
                    double N = initial_fine(a, i, j, k, test);
                    value += V * N;
                    if (test == 4)
                        momentum += V * N * sin(2 * PI * (i + .5) / a->fx);
                }
            n[slot * a->ns + k] = value / a->metric[p * a->ns + k].volume;
            m[slot * a->ns + k] = test == 4 ? momentum / a->metric[p * a->ns + k].volume
                                  : test    ? 2.5 * n[slot * a->ns + k]
                                            : 0;
        }
    }
}
/**
 * Integrate a conserved field with actual owner volumes and sum across MPI ranks.
 */
double a3_mass(A3 *a, const double *field) {
    double sum = 0;
    for (int p = 0; p < a->owned_count; p++)
        for (int k = 0; k < a->ns; k++)
            sum += field[p * a->ns + k] * a->metric[a->global[p] * a->ns + k].volume;
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}
/**
 * Collectively reject nonfinite/negative density, invalid vacuum momentum and out-of-range
 * marker states; never impose a density floor.
 */
int a3_valid(A3 *a, const double *n, const double *m) {
    int valid = 1;
    for (int t = 0; t < a->owned_count * a->ns; t++) {
        if (!isfinite(n[t]) || !isfinite(m[t]) || n[t] < 0 || (n[t] == 0 && m[t] != 0))
            valid = 0;
        if (a->marker && (!isfinite(a->marker[t]) || a->marker[t] < -2e-12 * a->density_scale ||
                          a->marker[t] > n[t] + 2e-12 * a->density_scale))
            valid = 0;
    }
    MPI_Allreduce(MPI_IN_PLACE, &valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return valid;
}
/**
 * Conservatively restrict/prolong density, parallel momentum and optional marker between old/new
 * active meshes, checking global transfer invariants.
 */
void a3_remap(A3 *old, A3 *next, double *n, double *m, double *new_n, double *new_m) {
    double mass = a3_mass(old, n), momentum = a3_mass(old, m), momentum_norm = 0;
    for (int p = 0; p < old->owned_count; p++)
        for (int k = 0; k < old->ns; k++)
            momentum_norm +=
                fabs(m[p * old->ns + k]) * old->metric[old->global[p] * old->ns + k].volume;
    MPI_Allreduce(MPI_IN_PLACE, &momentum_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    a3_prepare(old, n, m);
    next->material_min = old->material_min;
    next->material_max = old->material_max;
    for (int p = 0; p < next->owned_count; p++) {
        int id = next->global[p];
        for (int k = 0; k < next->ns; k++) {
            if (!a3_active(next, id, k))
                continue;
            A3Leaf l = a3_density_leaf(next, id, k);
            int s = a3_scale(l.level);
            double N = 0, M = 0;
            for (int v = 0; v < s; v++)
                for (int u = 0; u < s; u++) {
                    int i = s * l.i + u, j = s * l.j + v;
                    double V = next->g->cell[gid(next, 1, i, j, k)].volume;
                    N += V * a3_sample(old, 0, 1, i, j, k);
                    M += V * a3_sample(old, 1, 1, i, j, k);
                }
            new_n[p * next->ns + k] = N / next->metric[id * next->ns + k].volume;
            new_m[p * next->ns + k] = M / next->metric[id * next->ns + k].volume;
            if (next->marker) {
                double P = 0;
                for (int v = 0; v < s; v++)
                    for (int u = 0; u < s; u++) {
                        int i = s * l.i + u, j = s * l.j + v;
                        double V = next->g->cell[gid(next, 1, i, j, k)].volume;
                        P += V * a3_sample(old, 3, 1, i, j, k);
                    }
                next->marker[p * next->ns + k] = P / next->metric[id * next->ns + k].volume;
            }
        }
    }
    double err = fmax(fabs(a3_mass(next, new_n) - mass) / fmax(fabs(mass), 1e-300),
                      fabs(a3_mass(next, new_m) - momentum) / fmax(momentum_norm, 1e-300));
    next->remap_error = fmax(old->remap_error, err);
    next->boundary_total = old->boundary_total;
    if (!a3_valid(next, new_n, new_m) || err > 1e-11)
        a3_fail(next, "nonconservative or nonpositive volume remap");
}
/**
 * Release per-mesh arrays and PETSc objects; free immutable shared caches only when their
 * reference count reaches zero.
 */
void a3_free(A3 *a) {
    for (int f = 0; f < 3; f++)
        free(a->etd_weight[f]);
    free(a->macro_mask);
    free(a->macro_next);
    free(a->macro_field);
    free(a->macro_quiet);
    free(a->macro_born);
    free(a->macro_volume);
    free(a->macro_bg);
    for (int f = 0; f < 4; f++)
        free(a->macro_state[f]);
    a3_potential_free(a);
    for (int d = 0; d < 3; d++)
        free(a->potential_gradient[d]);
    free(a->parallel_face_geometry);
    free(a->face_blocks);
    free(a->field_requested);
    free(a->density_refined);
    free(a->tag_reason);
    free(a->seed_mask);
    free(a->quiet_since);
    free(a->born);
    free(a->tag_values);
    free(a->base_marker);
    free(a->refined);
    free(a->first);
    free(a->cover);
    free(a->slot);
    free(a->global);
    free(a->row);
    free(a->partition);
    free(a->base_velocity);
    free(a->velocity);
    free(a->log_n);
    free(a->child_n);
    free(a->child_m);
    free(a->child_velocity);
    free(a->virtual_offset);
    free(a->parallel_n);
    free(a->parallel_m);
    free(a->parallel_rate);
    free(a->rates);
    free(a->transverse_rates);
    free(a->leaf);
    free(a->metric);
    free(a->ready);
    for (int p = 0; p < a->peers; p++) {
        free(a->halo[p].recv);
        free(a->halo[p].send);
    }
    free(a->halo);
    free(a->requests);
    free(a->statuses);
    double *arrays[] = {a->base_n, a->base_m, a->sx_n, a->sy_n, a->amp_n};
    for (int i = 0; i < 5; i++)
        free(arrays[i]);
    if (a->cache_refs && --(*a->cache_refs) == 0) {
        free(a->cache_refs);
        free(a->fine_bg);
        free(a->response);
        free(a->root_volume);
        free(a->root_bg);
        free(a->offset_x);
        free(a->offset_y);
        for (int l = a->c.amr_levels == 3 ? -1 : 0; l < 2; l++) {
            free(a->parallel_length[l]);
            free(a->parallel_response[l]);
            for (int axis = 0; axis < 3; axis++)
                free(a->weights[l][axis]);
        }
    }
}
