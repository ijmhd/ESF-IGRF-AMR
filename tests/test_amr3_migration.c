#include "amr3_night_internal.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Abort the test communicator on a failed invariant; this stub supplies the migration module
 * test harness.
 */
_Noreturn void a3_fail(A3 *a, const char *s) {
    fprintf(stderr, "rank%d:%s\n", a->rank, s);
    MPI_Abort(MPI_COMM_WORLD, 2);
    abort();
}
/**
 * Allocate zeroed test storage or terminate the test communicator.
 */
void *a3_alloc(A3 *a, size_t n, size_t w) {
    void *p = calloc(n ? n : 1, w);
    if (!p)
        a3_fail(a, "allocation");
    return p;
}
/**
 * Generate a deterministic field/rank/index pattern so an incorrect ownership transfer is
 * observable.
 */
static double value(int field, size_t id) {
    return (field < 2 ? 1. : -1.) * ((field + 1) * 10000. + .25 * (id + 1));
}
/**
 * Compare migrated/retained values against the deterministic expected field pattern.
 */
static void verify(A3 *a, const int *target, int strict) {
    A3Night *q = a->night;
    size_t width = (size_t)2 * a->fx * a->ns;
    for (int f = 0; f < 4; f++)
        for (size_t x = 0; x < q->count; x++) {
            double *v = f < 2 ? q->ion[f] : q->momentum[f - 2];
            int old = x >= (size_t)a->row_first * width && x < (size_t)a->row_last * width;
            int next =
                x >= (size_t)target[a->rank] * width && x < (size_t)target[a->rank + 1] * width;
            if (old || next) {
                if (v[x] != value(f, x))
                    a3_fail(a, "changed or missing authoritative cell");
            } else if (strict && !isnan(v[x]))
                a3_fail(a, "unrelated cell overwritten");
        }
}
/**
 * Check simultaneous/cross-peer molecular migration, round trips, no-op byte identity and
 * preservation of old owned data.
 */
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    A3 a = {0};
    MPI_Comm_rank(MPI_COMM_WORLD, &a.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &a.size);
    a.fx = 8;
    a.ns = 6;
    a.ny = 8 * a.size;
    A3Night q = {0};
    a.night = &q;
    q.count = (size_t)a.fx * 2 * a.ny * a.ns;
    q.loss = .125;
    q.recomb = .25;
    a.partition = a3_alloc(&a, a.size + 1, sizeof(int));
    int *cuts = a3_alloc(&a, a.size + 1, sizeof(int)),
        *original = a3_alloc(&a, a.size + 1, sizeof(int));
    for (int r = 0; r <= a.size; r++) {
        original[r] = a.partition[r] = 8 * r;
        cuts[r] = 2 * (r + (3 * a.size * r * r) / (a.size * a.size));
    }
    a.row_first = a.partition[a.rank];
    a.row_last = a.partition[a.rank + 1];
    size_t width = (size_t)2 * a.fx * a.ns;
    for (int f = 0; f < 4; f++) {
        double **v = f < 2 ? q.ion + f : q.momentum + f - 2;
        *v = a3_alloc(&a, q.count, sizeof(double));
        for (size_t x = 0; x < q.count; x++)
            (*v)[x] = x >= (size_t)a.row_first * width && x < (size_t)a.row_last * width
                          ? value(f, x)
                          : NAN;
    }
    a3_night_repartition(&a, cuts);
    verify(&a, cuts, 1);
    memcpy(a.partition, cuts, (a.size + 1) * sizeof(int));
    a.row_first = cuts[a.rank];
    a.row_last = cuts[a.rank + 1];
    double *before = a3_alloc(&a, 4 * q.count, sizeof(double));
    for (int f = 0; f < 4; f++)
        memcpy(before + f * q.count, f < 2 ? q.ion[f] : q.momentum[f - 2],
               q.count * sizeof(double));
    a3_night_repartition(&a, cuts);
    for (int f = 0; f < 4; f++)
        if (memcmp(before + f * q.count, f < 2 ? q.ion[f] : q.momentum[f - 2],
                   q.count * sizeof(double)))
            a3_fail(&a, "no-op modified bytes");
    a3_night_repartition(&a, original);
    verify(&a, original, 0);
    memcpy(a.partition, original, (a.size + 1) * sizeof(int));
    a.row_first = original[a.rank];
    a.row_last = original[a.rank + 1];
    memcpy(cuts, original, (a.size + 1) * sizeof(int));
    if (a.size > 1)
        cuts[1] -= 2;
    a3_night_repartition(&a, cuts);
    verify(&a, cuts, 0);
    if (q.loss != .125 || q.recomb != .25)
        a3_fail(&a, "history modified");
    int ok = 1;
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!a.rank)
        printf(
            "MOLECULAR_OWNERSHIP_PASS ranks=%d simultaneous-cross-peer=1 round-trip=1 no-op-byte-identical=1 old-owned-preserved=1\n",
            a.size);
    free(before);
    for (int f = 0; f < 4; f++)
        free(f < 2 ? q.ion[f] : q.momentum[f - 2]);
    free(a.partition);
    free(cuts);
    free(original);
    MPI_Finalize();
    return 0;
}
