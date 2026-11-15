#ifndef ESF_AMR3_NIGHT_INTERNAL_H
#define ESF_AMR3_NIGHT_INTERNAL_H
#include "amr3.h"
#include "night.h"
#include "night_transport.h"
#include <stdint.h>
typedef struct A3Night {
    /* count is GLOBAL canonical cells, back_count is current OWNED cells.
     * ion/momentum retain full backing but only owned+refreshed ghosts are
     * valid between output gathers. Never read arbitrary remote backing. */
    size_t count, back_count;
    double *ion[2], *momentum[2], *backup[2], *backup_m[2];
    double *old_n[2], *old_m[2], *dn[2], *dm[2], *first_m[2], *rates[2];
    /* Frozen neutral-dependent transport/chemistry coefficients. electron
     * is recomputed from the current O+/molecular stage; it is not frozen. */
    double *coefficient[2], *nu[2], *rate, *electron;
    double *face[3];
    int face_first, face_last;
    double *transport_log[3]; /* Owned cells only; refreshed at each RHS stage. */
    double *etd_weight[2][3]; /* Same-attempt predictor/corrector coefficient reuse. */
    /* Owned cell averages at private levels -1,0,1, indexed by level+1.
     * They memoize the ORIGINAL volume-weighted sum/division, including on
     * the finest level. Validity is cleared whenever Ne is rebuilt. O+
     * pressure sampling is main-thread only, so these lazy writes do not
     * race with molecular worker kernels. No remote backing is cached. */
    double *electron_average[3];
    unsigned char *electron_ready[3];
    size_t electron_count[3];
    double alpha[2], loss, recomb, mass0, old_loss, old_recomb, max_charge_error, wall;
    double boundary_rate, boundary0, boundary_pending, boundary_total, transport_wall;
    uint64_t profile_hash;
} A3Night;
/* Global canonical-array offset of this rank's first owned molecular cell.
 * Each internal row corresponds to two canonical rows, with k contiguous. */
static inline size_t night_start(const A3 *a) {
    return (size_t)2 * a->row_first * a->fx * a->ns;
}
/* Number of owned canonical molecular cells, excluding ghost/backing-only
 * entries. Workspace arrays are indexed relative to night_start(a). */
static inline size_t night_owned(const A3 *a) {
    return (size_t)2 * (a->row_last - a->row_first) * a->fx * a->ns;
}
void a3_molecular_geometry(A3 *a);
void a3_molecular_workspace(A3 *a);
#endif
