/*
 * Legacy percentile detector and dispatch to the structure detector. The published teaching
 * cases explicitly select structure; legacy behavior remains available for comparison.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Same-time longitude reference removes common background evolution from tags.
 * Detection thresholds affect the mesh only, never density or conductivity.
 */
#include "amr3.h"
#include <math.h>
#include <stdlib.h>

/**
 * Compare two finite density values for the legacy longitude-percentile detector.
 */
static int order(const void *x, const void *y) {
    double a = *(const double *)x, b = *(const double *)y;
    return (a > b) - (a < b);
}
/**
 * Return the nonnegative periodic index i modulo n.
 */
static int wrap(int i, int n) {
    i %= n;
    return i < 0 ? i + n : i;
}
/**
 * Count a tag check and dispatch either the structure criterion or legacy percentile/depletion
 * criterion; both only request mesh changes.
 */
void a3_tags(A3 *a, double *n, double *m, unsigned char *mask, int *quiet) {
    a->tag_checks++;
    if (a->c.amr_structure) {
        a3_structure_tags(a, n, m, mask, quiet);
        return;
    }
    a3_prepare(a, n, m);
    double *ref = a3_alloc(a, (size_t)a->ny * a->ns, sizeof(double));
    double *fine_ref = a3_alloc(a, (size_t)a->fy * a->ns, sizeof(double));
    double *sorted = a3_alloc(a, a->fx, sizeof(double));
    unsigned char *core = a3_alloc(a, a->nb, 1), *want = a3_alloc(a, a->nb, 1);
    for (int j = 0; j < a->ny; j++)
        if (a->slot[a->first[j * a->nx]] >= 0)
            for (int k = 0; k < a->ns; k++) {
                for (int i = 0; i < a->nx; i++)
                    sorted[i] = a->base_n[(j * a->nx + i) * a->ns + k];
                qsort(sorted, a->nx, sizeof(double), order);
                ref[j * a->ns + k] = sorted[3 * (a->nx - 1) / 4];
            }
    /* A fine child's reference must be at the same fine apex label. Comparing
     * children at different heights with their parent's average mistakes the
     * unperturbed Chapman gradient for a new small-scale depletion. */
    for (int j = 0; j < a->fy; j++)
        if (a->c.amr_case == 1 || (j / 2 >= a->row_first - 1 && j / 2 <= a->row_last &&
                                   a->slot[a->first[(j / 2) * a->nx]] >= 0))
            for (int k = 0; k < a->ns; k++) {
                for (int i = 0; i < a->fx; i++)
                    sorted[i] = a3_sample(a, 0, 1, i, j, k);
                qsort(sorted, a->fx, sizeof(double), order);
                fine_ref[j * a->ns + k] = sorted[3 * (a->fx - 1) / 4];
            }
    for (int j = a->row_first; j < a->row_last; j++) {
        double peak = 0;
        for (int k = 0; k < a->ns; k++)
            peak = fmax(peak, ref[j * a->ns + k]);
        for (int i = 0; i < a->nx; i++) {
            int b = j * a->nx + i;
            double ratio = 1, gradient = 0;
            for (int k = 0; k < a->ns; k++) {
                double reference = ref[j * a->ns + k];
                TubeCell *g = a->g->cell + tube_cell_id(a->g, 2 * i, 2 * j, k);
                if (g->height < 200e3 || g->height > 1000e3 || reference <= .03 * peak)
                    continue;
                double lo = INFINITY, hi = 0;
                for (int c = 0; c < (a->refined[b] ? 4 : 1); c++) {
                    double denominator =
                        a->refined[b] ? fine_ref[(2 * j + (c >> 1)) * a->ns + k] : reference;
                    if (denominator <= .03 * peak)
                        continue;
                    double value = n[a->slot[a->first[b] + c] * a->ns + k] / denominator;
                    lo = fmin(lo, value);
                    hi = fmax(hi, value);
                }
                ratio = fmin(ratio, lo);
                gradient = fmax(gradient, hi - lo);
                for (int axis = 0; axis < 2; axis++) {
                    int il = i - (axis == 0), ir = i + (axis == 0), jl = j - (axis == 1),
                        jr = j + (axis == 1);
                    if (a->c.periodic_x) {
                        il = wrap(il, a->nx);
                        ir = wrap(ir, a->nx);
                    }
                    if (a->c.amr_case == 1) {
                        jl = wrap(jl, a->ny);
                        jr = wrap(jr, a->ny);
                    }
                    if (il < 0 || ir >= a->nx || jl < 0 || jr >= a->ny)
                        continue;
                    double rl = ref[jl * a->ns + k], rr = ref[jr * a->ns + k];
                    if (rl <= .03 * peak || rr <= .03 * peak)
                        continue;
                    double left = a->base_n[(jl * a->nx + il) * a->ns + k] / rl;
                    double right = a->base_n[(jr * a->nx + ir) * a->ns + k] / rr;
                    gradient = fmax(gradient, .5 * fabs(right - left));
                }
            }
            double dep = a->refined[b] ? a->c.amr_deplete_off : a->c.amr_deplete_on;
            double grad = a->refined[b] ? a->c.amr_grad_off : a->c.amr_grad_on;
            core[b] = (unsigned char)(ratio < dep || gradient > grad);
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, core, a->nb, MPI_UNSIGNED_CHAR, MPI_MAX, MPI_COMM_WORLD);
    for (int b = 0; b < a->nb; b++)
        if (core[b]) {
            int i = b % a->nx, j = b / a->nx;
            for (int v = -a->c.amr_buffer; v <= a->c.amr_buffer; v++)
                for (int u = -a->c.amr_buffer; u <= a->c.amr_buffer; u++) {
                    int x = i + u, y = j + v;
                    if (a->c.periodic_x)
                        x = wrap(x, a->nx);
                    if (a->c.amr_case == 1)
                        y = wrap(y, a->ny);
                    if (x >= 0 && x < a->nx && y >= 0 && y < a->ny)
                        want[y * a->nx + x] = 1;
                }
        }
    for (int b = 0; b < a->nb; b++) {
        if (want[b])
            quiet[b] = 0;
        else if (a->refined[b] && quiet[b] < a->c.amr_hold)
            quiet[b]++;
        else if (!a->refined[b])
            quiet[b] = 0;
        mask[b] = (unsigned char)(want[b] || (a->refined[b] && quiet[b] < a->c.amr_hold));
    }
    free(ref);
    free(fine_ref);
    free(sorted);
    free(core);
    free(want);
}
