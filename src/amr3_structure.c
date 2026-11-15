/*
 * Solution-based refinement detector, including represented-background residuals, velocity
 * curvature, seed protection and prospective-parent information loss. Detector thresholds never
 * modify plasma density or conductivity.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Local density DOFs, solution-based tags, and prospective-parent coarsening.
 * No percentile reference and no density/conductivity modification.
 */
#include "amr3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum {
    TAG_GRAD = 1,
    TAG_CURVE = 2,
    TAG_VELOCITY = 4,
    TAG_SEED = 8,
    TAG_LOSS = 16,
    TAG_BUFFER = 32,
    TAG_HOLD = 64,
    TAG_FIELD = 128,
    TAG_BACKGROUND = 256
};
/**
 * Return the unique evolved density backing owner for a column and k; account for both segmented
 * L1/L2 and macro L0 ownership.
 */
int a3_owner(const A3 *a, int p, int k) {
    int b = a->leaf[p].root;
    if (a->c.amr_levels == 3 && !a->macro_mask[a3_macro_index(a, b) * a->ns + k])
        return a3_macro_primary(a, b);
    return a->density_refined[b * a->ns + k] ? p : a->first[b];
}
/**
 * Test whether this backing slot is the actual evolved density owner at k.
 */
int a3_active(const A3 *a, int p, int k) {
    return a3_owner(a, p, k) == p;
}
/**
 * Return the physical density cell's level and labels, which can differ from the electric-
 * field/backing column level.
 */
A3Leaf a3_density_leaf(const A3 *a, int p, int k) {
    A3Leaf l = a->leaf[p];
    if (a->c.amr_levels == 3 && !a->macro_mask[a3_macro_index(a, l.root) * a->ns + k]) {
        l.level = -1;
        l.i = (l.root % a->nx) / 2;
        l.j = (l.root / a->nx) / 2;
        return l;
    }
    if (!a->density_refined[l.root * a->ns + k]) {
        l.level = 0;
        l.i = l.root % a->nx;
        l.j = l.root / a->nx;
    }
    return l;
}
/**
 * Select the segmented density mask or legacy column mask as required by the active case.
 */
const unsigned char *a3_density_mask(const A3 *a) {
    return a->c.amr_slab ? a->density_refined : a->refined;
}
/**
 * Apply boundary indexing on the detector's root/fine transverse stencil.
 */
static int index1(A3 *a, int axis, int q, int level) {
    int n = (axis ? a->ny : a->nx) << level;
    if (a->c.amr_case == 1 || (!axis && a->c.periodic_x)) {
        q %= n;
        return q < 0 ? q + n : q;
    }
    return q < 0 ? 0 : q >= n ? n - 1 : q;
}
/**
 * Evaluate or read a cached current E-cross-B drift averaged over canonical children at the
 * requested detector level.
 */
static V3 drift(A3 *a, int level, int i, int j, int k) {
    i = index1(a, 0, i, level);
    j = index1(a, 1, j, level);
    if (a->tag_drift_fine) {
        int first = a->tag_row_first << level;
        return (level ? a->tag_drift_fine
                      : a->tag_drift_root)[((j - first) * (a->nx << level) + i) * a->ns + k];
    }
    int s = 2 >> level;
    V3 out = {0};
    double V = 0;
    for (int v = 0; v < s; v++)
        for (int u = 0; u < s; u++) {
            int x = s * i + u, y = s * j + v;
            TubeCell *g = a->g->cell + tube_cell_id(a->g, x, y, k);
            double dx = a->fx * .5 *
                        (a->phi[y * (a->fx + 1) + x + 1] + a->phi[(y + 1) * (a->fx + 1) + x + 1] -
                         a->phi[y * (a->fx + 1) + x] - a->phi[(y + 1) * (a->fx + 1) + x]);
            double dy = a->fy * .5 *
                        (a->phi[(y + 1) * (a->fx + 1) + x] + a->phi[(y + 1) * (a->fx + 1) + x + 1] -
                         a->phi[y * (a->fx + 1) + x] - a->phi[y * (a->fx + 1) + x + 1]);
            V3 E = vadd(vmul(g->gx, -dx), vmul(g->gy, -dy));
            out = vadd(out, vmul(vcross(E, g->b), g->volume / g->B));
            V += g->volume;
        }
    return vmul(out, 1 / V);
}
/**
 * Estimate physical transverse cell spacing from contravariant metric magnitude, solely for
 * detector wavelength/velocity normalization.
 */
static double spacing(A3 *a, int level, int i, int j, int k, int axis) {
    int s = 2 >> level;
    TubeCell *g = a->g->cell + tube_cell_id(a->g, s * i, s * j, k);
    V3 v = axis ? g->gy : g->gx;
    return 1. / ((axis ? a->ny : a->nx) * (1 << level) * sqrt(vdot(v, v)));
}
/**
 * Compute normalized 3x3 density gradient/Hessian and drift-curvature indicators; optionally
 * remove the locally fitted represented-background component.
 */
static void indicator(A3 *a, int level, int i, int j, int k, double *gradient, double *curve,
                      double *speed, double *raw_gradient, double *raw_curve) {
    double N[3][3], norm = 1e-3 * a->density_scale;
    V3 U[3][3];
    for (int v = -1; v <= 1; v++)
        for (int u = -1; u <= 1; u++) {
            N[v + 1][u + 1] = a3_sample(a, 0, level, i + u, j + v, k);
            norm = fmax(norm, N[v + 1][u + 1]);
            if (a->c.amr_case != 1)
                U[v + 1][u + 1] = drift(a, level, i + u, j + v, k);
        }
    double dx = .5 * (N[1][2] - N[1][0]), dy = .5 * (N[2][1] - N[0][1]);
    double xx = N[1][2] - 2 * N[1][1] + N[1][0], yy = N[2][1] - 2 * N[1][1] + N[0][1];
    double xy = .25 * (N[2][2] - N[2][0] - N[0][2] + N[0][0]);
    *gradient = fmax(fabs(dx), fabs(dy)) / norm;
    *curve = sqrt(xx * xx + yy * yy + 2 * xy * xy) / norm;
    *speed = 0;
    *raw_gradient = *gradient;
    *raw_curve = *curve;
    if (a->c.amr_residual) {
        /* This is the background-shaped basis of the EXISTING conservative
         * prolongation, not an assumed time-dependent unseeded solution. */
        double B[3][3], R[3][3], num = 0, den = 0;
        for (int v = -1; v <= 1; v++)
            for (int u = -1; u <= 1; u++) {
                int x = index1(a, 0, i + u, level), y = index1(a, 1, j + v, level);
                B[v + 1][u + 1] = level ? a->fine_bg[tube_cell_id(a->g, x, y, k)]
                                        : a->root_bg[(y * a->nx + x) * a->ns + k];
            }
        const int X[] = {0, 2, 1, 1}, Y[] = {1, 1, 0, 2};
        for (int d = 0; d < 4; d++) {
            double bg = B[Y[d]][X[d]] - B[1][1];
            num += (N[Y[d]][X[d]] - N[1][1]) * bg;
            den += bg * bg;
        }
        double A = den > 0 ? num / den : 0;
        for (int v = 0; v < 3; v++)
            for (int u = 0; u < 3; u++)
                R[v][u] = (N[v][u] - N[1][1]) - A * (B[v][u] - B[1][1]);
        dx = .5 * (R[1][2] - R[1][0]);
        dy = .5 * (R[2][1] - R[0][1]);
        xx = R[1][2] + R[1][0];
        yy = R[2][1] + R[0][1];
        xy = .25 * (R[2][2] - R[2][0] - R[0][2] + R[0][0]);
        *gradient = fmax(fabs(dx), fabs(dy)) / norm;
        *curve = sqrt(xx * xx + yy * yy + 2 * xy * xy) / norm;
    }
    if (a->c.amr_case != 1) {
        V3 vx = vadd(vsub(U[1][2], vmul(U[1][1], 2)), U[1][0]);
        V3 vy = vadd(vsub(U[2][1], vmul(U[1][1], 2)), U[0][1]);
        V3 vm = vmul(vsub(vsub(U[2][2], U[2][0]), vsub(U[0][2], U[0][0])), .25);
        double h = fmin(spacing(a, level, i, j, k, 0), spacing(a, level, i, j, k, 1));
        *speed = a->c.amr_interval * sqrt(vdot(vx, vx) + vdot(vy, vy) + 2 * vdot(vm, vm)) / h;
    }
}
/**
 * Compute/cache the root-cell field-line-integrated conductivity tensor used to detect elliptic-
 * coefficient curvature.
 */
static void integrated(A3 *a, int i, int j, double out[4]) {
    i = index1(a, 0, i, 0);
    j = index1(a, 1, j, 0);
    if (a->tag_tensor) {
        memcpy(out, a->tag_tensor + 4 * ((j - a->tag_row_first) * a->nx + i), 4 * sizeof(double));
        return;
    }
    memset(out, 0, 4 * sizeof(double));
    for (int k = 0; k < a->ns; k++)
        for (int c = 0; c < 4; c++) {
            int x = 2 * i + (c & 1), y = 2 * j + (c >> 1), id = tube_cell_id(a->g, x, y, k);
            double N = a3_sample(a, 0, 1, x, y, k), *r = a->response + 12 * id;
            for (int d = 0; d < 4; d++)
                out[d] += (r[d] * N + r[4 + d]);
        }
}
/**
 * Measure normalized transverse second differences of the complete integrated tensor, including
 * independent Hall cross terms.
 */
static double coefficient_curve(A3 *a, int i, int j) {
    double C[4], L[4], R[4], D[4], U[4];
    integrated(a, i, j, C);
    integrated(a, i - 1, j, L);
    integrated(a, i + 1, j, R);
    integrated(a, i, j - 1, D);
    integrated(a, i, j + 1, U);
    double num = 0, den = 0;
    for (int d = 0; d < 4; d++) {
        double x = L[d] - 2 * C[d] + R[d], y = D[d] - 2 * C[d] + U[d];
        num += x * x + y * y;
        den +=
            fmax(C[d] * C[d], fmax(fmax(L[d] * L[d], R[d] * R[d]), fmax(D[d] * D[d], U[d] * U[d])));
    }
    return den > 0 ? sqrt(num / den) : 0;
}
/**
 * Check declared seed/envelope sampling requirements and count locations that remain unresolved
 * even on the finest available mesh.
 */
static int mode_needs_fine(A3 *a, int i, int j, int k, double lambda, int envelope) {
    if (!(lambda > 0))
        return 0;
    double hx = spacing(a, 0, i, j, k, 0), hy = spacing(a, 0, i, j, k, 1);
    int need = hx > lambda / 32;
    if (envelope)
        need |= hy > a->c.seed_width / 8;
    if (.5 * hx > lambda / 32 * (1 + 1e-10) ||
        (envelope && .5 * hy > a->c.seed_width / 8 * (1 + 1e-10)))
        a->unresolved_seed++;
    return need;
}
/**
 * Collectively form structure-based refinement requests, prospective-parent loss tests, seed
 * guard, slab promotion, buffers and hysteresis. No unperturbed control solution is evolved by
 * this detector.
 */
void a3_structure_tags(A3 *a, double *n, double *m, unsigned char *mask, int *quiet) {
    double start = wall3();
    a3_prepare(a, n, m);
    if (a->c.amr_case != 1)
        a3_potential(a);
    if (a->c.amr_case != 1) {
        a->tag_row_first = a->row_first ? a->row_first - 1 : 0;
        a->tag_row_last = a->row_last < a->ny ? a->row_last + 1 : a->ny;
        int rows = a->tag_row_last - a->tag_row_first;
        V3 *fine = a3_alloc(a, (size_t)4 * rows * a->nx * a->ns, sizeof(V3));
        V3 *root = a3_alloc(a, (size_t)rows * a->nx * a->ns, sizeof(V3));
        double *tensor = a3_alloc(a, (size_t)4 * rows * a->nx, sizeof(double));
        for (int j = 2 * a->tag_row_first; j < 2 * a->tag_row_last; j++)
            for (int i = 0; i < a->fx; i++)
                for (int k = 0; k < a->ns; k++)
                    fine[((j - 2 * a->tag_row_first) * a->fx + i) * a->ns + k] =
                        drift(a, 1, i, j, k);
        for (int j = a->tag_row_first; j < a->tag_row_last; j++)
            for (int i = 0; i < a->nx; i++) {
                integrated(a, i, j, tensor + 4 * ((j - a->tag_row_first) * a->nx + i));
                for (int k = 0; k < a->ns; k++) {
                    V3 sum = {0};
                    double volume = 0;
                    for (int v = 0; v < 2; v++)
                        for (int u = 0; u < 2; u++) {
                            int x = 2 * i + u, y = 2 * j + v;
                            double V = a->g->cell[tube_cell_id(a->g, x, y, k)].volume;
                            sum = vadd(
                                sum,
                                vmul(fine[((y - 2 * a->tag_row_first) * a->fx + x) * a->ns + k],
                                     V));
                            volume += V;
                        }
                    root[((j - a->tag_row_first) * a->nx + i) * a->ns + k] = vmul(sum, 1 / volume);
                }
            }
        a->tag_drift_fine = fine;
        a->tag_drift_root = root;
        a->tag_tensor = tensor;
    }
    size_t cells = (size_t)a->nb * a->ns;
    unsigned char *core = a3_alloc(a, cells, 1), *want = a3_alloc(a, cells, 1);
    memset(a->tag_reason, 0, cells * sizeof(unsigned short));
    memset(a->tag_values, 0, cells * 6 * sizeof(double));
    a->unresolved_seed = 0;
    memset(a->field_requested, 0, a->nb);
    for (int j = a->row_first; j < a->row_last; j++)
        for (int i = 0; i < a->nx; i++) {
            int b = j * a->nx + i;
            double coeff = a->c.amr_case == 1 ? 0 : coefficient_curve(a, i, j);
            a->field_requested[b] =
                (unsigned char)(coeff > (a->refined[b] ? .5 : 1) * a->c.amr_sigma_on);
            for (int k = 0; k < a->ns; k++) {
                int t = b * a->ns + k, old = a->density_refined[t];
                double G, C, V, rawG, rawC;
                indicator(a, 0, i, j, k, &G, &C, &V, &rawG, &rawC);
                unsigned short reason = 0;
                if (G > (old ? a->c.amr_grad_off : a->c.amr_grad_on))
                    reason |= TAG_GRAD;
                if (C > (old ? a->c.amr_curve_off : a->c.amr_curve_on))
                    reason |= TAG_CURVE;
                if (V > (old ? .5 : 1) * a->c.amr_velocity_on)
                    reason |= TAG_VELOCITY;
                double loss = 0;
                if (old) {
                    /* True fine structure must survive restriction. In particular,
                     * canceling child details cannot disappear in the parent mean. */
                    double norm = 1e-3 * a->density_scale;
                    for (int c = 0; c < 4; c++)
                        norm = fmax(norm, a3_sample(a, 0, 1, 2 * i + (c & 1), 2 * j + (c >> 1), k));
                    for (int c = 0; c < 4; c++) {
                        double q = a3_sample(a, 0, 1, 2 * i + (c & 1), 2 * j + (c >> 1), k);
                        double fg, fc, fv, rg, rc;
                        indicator(a, 1, 2 * i + (c & 1), 2 * j + (c >> 1), k, &fg, &fc, &fv, &rg,
                                  &rc);
                        if (fg > a->c.amr_grad_on)
                            reason |= TAG_GRAD;
                        if (fc > a->c.amr_curve_on)
                            reason |= TAG_CURVE;
                        if (fv > a->c.amr_velocity_on)
                            reason |= TAG_VELOCITY;
                        /* Same background-preserving parent prolongation as the
                         * production remap, exposed independently of active leaves. */
                        double predicted = a3_proposed_child(a, b, c, k, 0);
                        loss = fmax(loss, fabs(q - predicted) / norm);
                    }
                    double allowed = .005;
                    if (a->marker && a->base_n[t] > 0 && a->base_marker[t] / a->base_n[t] >= 1e-3)
                        allowed = fmin(allowed, .1 * a->c.amr_seed_amplitude);
                    if (loss > allowed)
                        reason |= TAG_LOSS;
                }
                if (a->marker && a->base_n[t] > 0 && a->base_marker[t] / a->base_n[t] >= 1e-3 &&
                    mode_needs_fine(a, i, j, k, a->seed_lambda, 1))
                    reason |= TAG_SEED;
                if ((reason & (TAG_GRAD | TAG_CURVE)) &&
                    mode_needs_fine(a, i, j, k, a->branch_lambda, 0))
                    reason |= TAG_SEED;
                if (a->field_requested[b])
                    reason |= TAG_FIELD;
                if (a->c.amr_residual && !(reason & (TAG_GRAD | TAG_CURVE)) &&
                    (rawG > (old ? a->c.amr_grad_off : a->c.amr_grad_on) ||
                     rawC > (old ? a->c.amr_curve_off : a->c.amr_curve_on)))
                    reason |= TAG_BACKGROUND;
                a->tag_reason[t] = reason;
                core[t] = (unsigned char)((reason & ~(TAG_FIELD | TAG_BACKGROUND)) != 0);
                double *value = a->tag_values + 6 * (size_t)t;
                value[0] = G;
                value[1] = C;
                value[2] = V;
                value[3] = loss;
                value[4] = rawG;
                value[5] = rawC;
            }
        }
    MPI_Allreduce(MPI_IN_PLACE, core, (int)cells, MPI_UNSIGNED_CHAR, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, a->tag_reason, (int)cells, MPI_UNSIGNED_SHORT, MPI_BOR,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, a->tag_values, (int)(6 * cells), MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, a->field_requested, a->nb, MPI_UNSIGNED_CHAR, MPI_MAX,
                  MPI_COMM_WORLD);
    unsigned long unresolved = a->unresolved_seed;
    MPI_Allreduce(&unresolved, &a->unresolved_seed, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    /* Two fine-cell template supports = one root cell. Additional displacement
     * is bounded by an earlier regrid at half a root cell. No all-height OR. */
    int buffer = a->c.amr_buffer;
    for (int j = 0; j < a->ny; j++)
        for (int i = 0; i < a->nx; i++)
            for (int k = 0; k < a->ns; k++)
                if (core[(j * a->nx + i) * a->ns + k])
                    for (int v = -buffer; v <= buffer; v++)
                        for (int u = -buffer; u <= buffer; u++) {
                            int x = i + u, y = j + v;
                            if (a->c.periodic_x || a->c.amr_case == 1)
                                x = index1(a, 0, x, 0);
                            if (a->c.amr_case == 1)
                                y = index1(a, 1, y, 0);
                            if (x >= 0 && x < a->nx && y >= 0 && y < a->ny) {
                                int t = (y * a->nx + x) * a->ns + k;
                                want[t] = 1;
                                if (!core[t])
                                    a->tag_reason[t] |= TAG_BUFFER;
                            }
                        }
    int slab = a->c.amr_slab ? a->c.amr_slab : a->ns;
    for (int b = 0; b < a->nb; b++)
        for (int k0 = 0; k0 < a->ns; k0 += slab) {
            int requested = 0;
            for (int k = k0; k < k0 + slab; k++)
                requested |= want[b * a->ns + k];
            for (int k = k0; k < k0 + slab; k++) {
                int cell = b * a->ns + k, index = a->c.amr_slab ? cell : b,
                    old = a->density_refined[cell];
                if (requested) {
                    a->quiet_since[index] = -1;
                    quiet[index] = 0;
                    mask[index] = 1;
                } else if (old) {
                    if (a->quiet_since[index] < 0)
                        a->quiet_since[index] = a->tag_time;
                    int wait = a->tag_time - a->quiet_since[index] < a->c.amr_quiet_s ||
                               a->tag_time - a->born[index] < a->c.amr_quiet_s;
                    mask[index] = (unsigned char)wait;
                    if (wait)
                        a->tag_reason[cell] |= TAG_HOLD;
                } else {
                    mask[index] = 0;
                    a->quiet_since[index] = -1;
                }
            }
        }
    if (a->c.amr_levels == 3)
        a3_macro_tags(a, mask);
    a->tag_elapsed += wall3() - start;
    free(core);
    free(want);
    free(a->tag_drift_fine);
    free(a->tag_drift_root);
    free(a->tag_tensor);
    a->tag_drift_fine = a->tag_drift_root = NULL;
    a->tag_tensor = NULL;
}
/**
 * Open one detector diagnostic output or abort on failure.
 */
static FILE *out(A3 *a, const char *name, const char *mode) {
    char p[1200];
    snprintf(p, sizeof(p), "%s/%s", a->c.out, name);
    FILE *f = fopen(p, mode);
    if (!f)
        a3_fail(a, "structure audit open");
    return f;
}
/**
 * Write actual density masks, tag reason bits/values, marker diagnostics and frame-to-mesh
 * bindings for plotting.
 */
void a3_structure_output(A3 *a, const double *n, const double *marker, int frame, double time) {
    if (!a->c.amr_structure || a->rank)
        return;
    size_t cells = (size_t)a->nb * a->ns;
    char name[100];
    FILE *f;
    snprintf(name, sizeof(name), "density_mask_%04d.bin", a->mesh_epoch);
    f = out(a, name, "wb");
    if (fwrite(a->density_refined, 1, cells, f) != cells)
        a3_fail(a, "density mask write");
    fclose(f);
    snprintf(name, sizeof(name), "tag_reasons_%04d.bin", frame);
    f = out(a, name, "wb");
    if (fwrite(a->tag_reason, sizeof(unsigned short), cells, f) != cells)
        a3_fail(a, "tag reasons write");
    fclose(f);
    snprintf(name, sizeof(name), "tag_values_%04d.bin", frame);
    f = out(a, name, "wb");
    if (fwrite(a->tag_values, sizeof(double), 6 * cells, f) != 6 * cells)
        a3_fail(a, "tag values write");
    fclose(f);
    if (marker) {
        snprintf(name, sizeof(name), "marker_%04d.bin", frame);
        f = out(a, name, "wb");
        if (fwrite(marker, sizeof(double), (size_t)a->nleaf * a->ns, f) != (size_t)a->nleaf * a->ns)
            a3_fail(a, "marker write");
        fclose(f);
    }
    if (!frame) {
        f = out(a, "structure_frames.csv", "w");
        fputs(
            "frame,time_s,tag_time_s,mesh_epoch,active_density_cells,backing_slots,field_unknowns,density_mask,tag_reasons,tag_values,marker\n",
            f);
        fclose(f);
        f = out(a, "structure_metadata.json", "w");
        fprintf(
            f,
            "{\"schema\":3,\"density_indicator\":\"%s\",\"root_shape\":[%d,%d,%d],\"slab_cells\":%d,\"gradient_on\":%.17g,\"gradient_off\":%.17g,\"curvature_on\":%.17g,\"curvature_off\":%.17g,\"velocity_on\":%.17g,\"coefficient_on\":%.17g,\"quiet_seconds\":%.17g,\"density_reference_scale\":%.17g,\"seed_wavelength_m\":%.17g,\"branch_wavelength_m\":%.17g,\"seed_amplitude_budget\":%.17g,\"storage\":\"union-column backing; zero-volume inactive slots; only active FV unknowns evolved; storage is not yet compacted\",\"tag_value_order\":[\"parent_gradient\",\"parent_curvature\",\"parent_velocity_indicator\",\"child_reconstruction_loss\",\"raw_parent_gradient\",\"raw_parent_curvature\"],\"tag_bits\":{\"density_gradient\":1,\"density_curvature\":2,\"velocity_curvature\":4,\"mode_resolution\":8,\"coarsening_loss\":16,\"motion_buffer\":32,\"coarsening_hold\":64,\"field_coefficient\":128,\"represented_background\":256}}\n",
            a->c.amr_residual ? "represented" : "total", a->nx, a->ny, a->ns, a->c.amr_slab,
            a->c.amr_grad_on, a->c.amr_grad_off, a->c.amr_curve_on, a->c.amr_curve_off,
            a->c.amr_velocity_on, a->c.amr_sigma_on, a->c.amr_quiet_s, a->density_scale,
            a->seed_lambda, a->branch_lambda, a->c.amr_seed_amplitude);
        fclose(f);
    }
    f = out(a, "structure_frames.csv", "a");
    fprintf(
        f,
        "%d,%.17g,%.17g,%d,%d,%d,%d,density_mask_%04d.bin,tag_reasons_%04d.bin,tag_values_%04d.bin,",
        frame, time, a->tag_time, a->mesh_epoch, a->active_cells, a->nleaf * a->ns, a->ndof,
        a->mesh_epoch, frame, frame);
    if (marker)
        fprintf(f, "marker_%04d.bin", frame);
    fputc('\n', f);
    fclose(f);
    double volume = 0, fine_volume = 0, marker_mass = 0, low = 1, high = 0;
    unsigned long counts[9] = {0};
    for (size_t q = 0; q < cells; q++) {
        volume += a->root_volume[q];
        if (a->density_refined[q])
            fine_volume += a->root_volume[q];
        for (int bit = 0; bit < 9; bit++)
            counts[bit] += (a->tag_reason[q] >> bit) & 1;
    }
    if (marker)
        for (int p = 0; p < a->nleaf; p++)
            for (int k = 0; k < a->ns; k++)
                if (a3_active(a, p, k)) {
                    int q = p * a->ns + k;
                    marker_mass += marker[q] * a->metric[q].volume;
                    double chi = n[q] > 0 ? marker[q] / n[q] : 0;
                    low = fmin(low, chi);
                    high = fmax(high, chi);
                }
    f = out(a, "structure_status.json", "w");
    fprintf(
        f,
        "{\"time_s\":%.17g,\"last_tag_time_s\":%.17g,\"active_density_cells\":%d,\"backing_slots\":%d,\"uniform_finest_cells\":%d,\"refined_physical_volume_fraction\":%.17g,\"tag_seconds\":%.9g,\"unresolved_mode_locations\":%lu,\"marker_mass\":%.17g,\"marker_fraction_min\":%.17g,\"marker_fraction_max\":%.17g,\"reason_counts\":[%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu],\"refined_roots_per_along_index\":[",
        time, a->tag_time, a->active_cells, a->nleaf * a->ns, 4 * a->nb * a->ns,
        fine_volume / volume, a->tag_elapsed, a->unresolved_seed, marker_mass, low, high, counts[0],
        counts[1], counts[2], counts[3], counts[4], counts[5], counts[6], counts[7], counts[8]);
    for (int k = 0; k < a->ns; k++) {
        int count = 0;
        for (int b = 0; b < a->nb; b++)
            count += a->density_refined[b * a->ns + k];
        fprintf(f, "%s%d", k ? "," : "", count);
    }
    fputs("]}\n", f);
    fclose(f);
}
