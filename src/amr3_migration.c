/*
 * Molecular old/new slab ownership intersections. Separate receive buffers preserve every old
 * owned value until O+ remapping and all sends finish.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* SPDX-License-Identifier: GPL-3.0-only
 * Redistribute only molecular cells whose owner changes. Old owned values stay
 * valid until O+ remap is finished; separate receives prevent send overwrite.
 */
#include "amr3_night_internal.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>

/**
 * Move only molecular cells whose owner changes using old/new slab intersections and separate
 * receive storage. Identical intervals are a byte-preserving no-op; old owned values survive
 * until O+ remap completes.
 */
void a3_night_repartition(A3 *a, const int *cuts) {
    A3Night *q = a->night;
    if (!q || !memcmp(a->partition, cuts, (a->size + 1) * sizeof(int)))
        return;
    if (cuts[0] != 0 || cuts[a->size] != a->ny)
        a3_fail(a, "molecular redistribution endpoints");
    for (int r = 0; r < a->size; r++)
        if (cuts[r] >= cuts[r + 1])
            a3_fail(a, "molecular redistribution empty owner");
    if (cuts[a->rank] == a->row_first && cuts[a->rank + 1] == a->row_last)
        return;
    size_t width = (size_t)2 * a->fx * a->ns;
    size_t first = (size_t)cuts[a->rank] * width;
    size_t count = (size_t)(cuts[a->rank + 1] - cuts[a->rank]) * width;
    if (count > INT_MAX)
        a3_fail(a, "molecular redistribution MPI count overflow");
    double *received = a3_alloc(a, 4 * count, sizeof(double));
    MPI_Request *requests = a3_alloc(a, 8 * a->size, sizeof(MPI_Request));
    MPI_Status *statuses = a3_alloc(a, 8 * a->size, sizeof(MPI_Status));
    int nr = 0;
    for (int field = 0; field < 4; field++) {
        double *v = field < 2 ? q->ion[field] : q->momentum[field - 2];
        for (int peer = 0; peer < a->size; peer++) {
            int lo = cuts[a->rank] > a->partition[peer] ? cuts[a->rank] : a->partition[peer];
            int hi = cuts[a->rank + 1] < a->partition[peer + 1] ? cuts[a->rank + 1]
                                                                : a->partition[peer + 1];
            if (hi > lo) {
                size_t n = (size_t)(hi - lo) * width, offset = (size_t)lo * width;
                if (n > INT_MAX)
                    a3_fail(a, "molecular receive MPI count overflow");
                if (peer == a->rank)
                    memcpy(received + field * count + (offset - first), v + offset,
                           n * sizeof(double));
                else
                    MPI_Irecv(received + field * count + (offset - first), (int)n, MPI_DOUBLE, peer,
                              190 + field, MPI_COMM_WORLD, requests + nr++);
            }
            if (peer == a->rank)
                continue;
            lo = a->row_first > cuts[peer] ? a->row_first : cuts[peer];
            hi = a->row_last < cuts[peer + 1] ? a->row_last : cuts[peer + 1];
            if (hi > lo) {
                size_t n = (size_t)(hi - lo) * width;
                if (n > INT_MAX)
                    a3_fail(a, "molecular send MPI count overflow");
                MPI_Isend(v + (size_t)lo * width, (int)n, MPI_DOUBLE, peer, 190 + field,
                          MPI_COMM_WORLD, requests + nr++);
            }
        }
    }
    MPI_Waitall(nr, requests, statuses);
    for (int field = 0; field < 4; field++) {
        double *v = field < 2 ? q->ion[field] : q->momentum[field - 2];
        memcpy(v + first, received + field * count, count * sizeof(double));
    }
    free(statuses);
    free(requests);
    free(received);
}
