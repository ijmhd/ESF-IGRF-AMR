/*
 * Uniform-grid reference solver. Density is normalized by 1e12 m^-3 and momentum is
 * n*u_parallel. This path is retained for static-grid and backend comparisons.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define PI 3.14159265358979323846
/**
 * Load IGRF, degree-one dipole or the explicit aligned-dipole fixture from the selected magnetic
 * policy.
 */
int solver3_load_field(const Config3 *c, Igrf *field) {
    int code = igrf_load(field, c->coeff, c->year, c->degree);
    if (!code && c->aligned_dipole) {
        memset(field->g, 0, sizeof(field->g));
        memset(field->h, 0, sizeof(field->h));
        field->degree = 1;
        field->g[1][0] = -c->dipole_equator_nT;
    }
    return code;
}
/**
 * Evaluate the selected idealized single/double Chapman O+ altitude profile in normalized units.
 */
double solver3_initial_density(const Config3 *c, double h) {
    /* TELF Eq.13: both Chapman populations are transported O+, not a fixed
     * E-region species. N is normalized to 1e12 m^-3. */
    return c->background_profile ? chapman(h, 438e3, 80e3) + .01 * chapman(h, 600e3, 200e3)
                                 : chapman(h, 438e3, 70e3);
}
/**
 * Report a fatal rank-local model error and abort MPI when enabled.
 */
_Noreturn void fail3(Solver3 *s, const char *message) {
    fprintf(stderr, "rank %d: %s\n", s->rank, message);
#ifdef ESF_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    exit(1);
}
/**
 * Allocate zero-initialized reference-solver storage and fail on allocation errors.
 */
void *alloc3(Solver3 *s, size_t count, size_t width) {
    void *p = calloc(count, width);
    if (!p)
        fail3(s, "allocation failed");
    return p;
}
/**
 * Return a monotonic wall clock for timing; this is independent of model physical time.
 */
double wall3(void) {
    struct timespec t;
#ifdef __APPLE__
    clock_gettime(CLOCK_UPTIME_RAW, &t);
#else
    clock_gettime(CLOCK_MONOTONIC, &t);
#endif
    return t.tv_sec + 1e-9 * t.tv_nsec;
}
/**
 * Collect distributed uniform owned state slabs into replicated global arrays when needed.
 */
void solver3_gather(Solver3 *s, double *x) {
#ifdef ESF_MPI
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DOUBLE, x, s->counts, s->offsets, MPI_DOUBLE,
                   MPI_COMM_WORLD);
#else
    (void)s;
    (void)x;
#endif
}
/* Exchange only the two-row MUSCL/PDM4 dependency halo, with both fields in
 * one message per peer. Full root gathers are reserved for output. */
/**
 * Exchange uniform-grid density/momentum halos while retaining owned and ghost consistency.
 */
static void exchange_state(Solver3 *s, double *n, double *m) {
#ifdef ESF_MPI
    if (s->size == 1)
        return;
    MPI_Request *requests = s->halo_requests;
    int q = 0, ns = s->c.grid.ns;
    for (int k = 0; k < s->halo_peers; k++) {
        HaloPeer3 *h = s->halo + k;
        if (h->recv_count)
            MPI_Irecv(h->recv, 2 * h->recv_count * ns, MPI_DOUBLE, h->peer, 73, MPI_COMM_WORLD,
                      requests + q++);
    }
    for (int k = 0; k < s->halo_peers; k++) {
        HaloPeer3 *h = s->halo + k;
        int len = h->send_count * ns, off = h->send_first * ns;
        if (len) {
            memcpy(h->send, n + off, len * sizeof(double));
            memcpy(h->send + len, m + off, len * sizeof(double));
            MPI_Isend(h->send, 2 * len, MPI_DOUBLE, h->peer, 73, MPI_COMM_WORLD, requests + q++);
        }
    }
    MPI_Waitall(q, requests, s->halo_statuses);
    for (int k = 0; k < s->halo_peers; k++) {
        HaloPeer3 *h = s->halo + k;
        int len = h->recv_count * ns, off = h->recv_first * ns;
        if (len) {
            memcpy(n + off, h->recv, len * sizeof(double));
            memcpy(m + off, h->recv + len, len * sizeof(double));
        }
    }
#else
    (void)s;
    (void)n;
    (void)m;
#endif
}
/**
 * Broadcast a potentially large raw array in MPI-count-safe chunks.
 */
static void broadcast(Solver3 *s, void *x, size_t bytes) {
#ifdef ESF_MPI
    if (bytes > 2147483647)
        fail3(s, "geometry exceeds the current MPI buffer size");
    MPI_Bcast(x, (int)bytes, MPI_BYTE, 0, MPI_COMM_WORLD);
#else
    (void)s;
    (void)x;
    (void)bytes;
#endif
}
/**
 * Construct uniform geometry, prescribed background, decomposition, reconstruction caches and
 * potential backend.
 */
void solver3_create(Solver3 *s, Config3 c, int rank, int size) {
    memset(s, 0, sizeof(*s));
    s->c = c;
    s->rank = rank;
    s->size = size;
    s->wall0 = wall3();
    if (solver3_load_field(&c, &s->field))
        fail3(s, "invalid field coefficients/date");
    if (tube_apex(&s->field, geodetic_position(22.3, 114.2, 350e3), &s->target_apex))
        fail3(s, "Hong Kong field line did not reach an apex");
    if (!rank) {
        fprintf(stdout, "Tracing %d x %d complete conjugate flux tubes...\n", c.grid.nx + 1,
                c.grid.ny + 1);
        fflush(stdout);
        int code = tube_geometry(&s->g, &s->field, c.grid);
        if (code) {
            fprintf(stderr, "geometry error %d\n", code);
            fail3(s, "invalid field-aligned geometry");
        }
        tube_conductivity(&s->g, c.e_peak);
        if (c.collision_profile[0] && tube_conductivity_table(&s->g, c.e_peak, c.collision_profile))
            fail3(
                s,
                "invalid collision table: require ordered height_m,nu_O_s1,nu_e_s1 with positive frequencies and full altitude coverage");
    }
    if (rank) {
        s->g.c = c.grid;
        s->g.node =
            alloc3(s, (size_t)(c.grid.nx + 1) * (c.grid.ny + 1) * (c.grid.ns + 1), sizeof(V3));
        s->g.cell = alloc3(s, (size_t)c.grid.nx * c.grid.ny * c.grid.ns, sizeof(TubeCell));
        for (int a = 0; a < 3; a++)
            s->g.face[a] = alloc3(s, tube_face_count(&s->g, a), sizeof(TubeFace));
    }
    broadcast(s, s->g.node,
              (size_t)(c.grid.nx + 1) * (c.grid.ny + 1) * (c.grid.ns + 1) * sizeof(V3));
    broadcast(s, s->g.cell, (size_t)c.grid.nx * c.grid.ny * c.grid.ns * sizeof(TubeCell));
    for (int a = 0; a < 3; a++)
        broadcast(s, s->g.face[a], tube_face_count(&s->g, a) * sizeof(TubeFace));
    broadcast(s, &s->g.max_alignment_error, sizeof(double));
    broadcast(s, &s->g.max_parallel_electric_metric_error, sizeof(double));
    broadcast(s, &s->g.min_volume, sizeof(double));
    broadcast(s, &s->g.max_magnetic_flux_jump, sizeof(double));
    if (c.exb_pdm4) {
        for (int a = 0; a < 3; a++)
            s->pdm_weights[a] = alloc3(s, 4 * (size_t)tube_face_count(&s->g, a), sizeof(double));
        if (tube_pdm_weights(&s->g, c.periodic_x, s->pdm_weights))
            fail3(s, "PDM4 physical-volume moment construction failed");
    }
    if (c.background_east != 0 || c.seed_east > 0) {
        /* Use a fixed physical reference, independent of grid centers. A
         * nearest-cell calibration changes both forcing and seed width during
         * a refinement study. Derive the perpendicular dual basis directly on
         * the smooth apex surface at the exact seed height and central label. */
        V3 ref, xl, xr, yl, yr;
        double dl = 1e-4, dh = 10;
        if (tube_apex_label(&s->field, c.grid.lon + c.seed_offset_deg, c.seed_height, &ref) ||
            tube_apex_label(&s->field, c.grid.lon + c.seed_offset_deg - dl, c.seed_height, &xl) ||
            tube_apex_label(&s->field, c.grid.lon + c.seed_offset_deg + dl, c.seed_height, &xr) ||
            tube_apex_label(&s->field, c.grid.lon + c.seed_offset_deg, c.seed_height - dh, &yl) ||
            tube_apex_label(&s->field, c.grid.lon + c.seed_offset_deg, c.seed_height + dh, &yr))
            fail3(s, "cannot calibrate fixed apex reference");
        V3 ex = vmul(vsub(xr, xl), c.grid.span / (2 * dl)),
           ey = vmul(vsub(yr, yl), (c.grid.hmax - c.grid.hmin) / (2 * dh)),
           b = vunit(igrf_field(&s->field, ref));
        V3 gx = vmul(vcross(ey, b), 1 / vdot(ex, vcross(ey, b)));
        double lon = atan2(ref.y, ref.x), projection = vdot(gx, (V3){-sin(lon), cos(lon), 0});
        if (projection <= 0)
            fail3(s, "invalid reference direction for imposed eastward field");
        s->east_scale_m = 1 / projection;
        s->background_potential_gradient = -c.background_east * s->east_scale_m;
    }
    int columns = c.grid.nx * c.grid.ny;
    s->first = rank * (columns / size) + (rank < columns % size ? rank : columns % size);
    s->count = columns / size + (rank < columns % size);
    if (!s->count)
        fail3(s, "more MPI ranks than flux tubes");
    s->n = s->count * c.grid.ns;
    s->total = columns * c.grid.ns;
    s->counts = alloc3(s, size, sizeof(int));
    s->offsets = alloc3(s, size, sizeof(int));
    for (int r = 0; r < size; r++) {
        s->counts[r] = (columns / size + (r < columns % size)) * c.grid.ns;
        s->offsets[r] =
            (r * (columns / size) + (r < columns % size ? r : columns % size)) * c.grid.ns;
    }
    s->coefficient_counts = alloc3(s, size, sizeof(int));
    s->coefficient_offsets = alloc3(s, size, sizeof(int));
    s->halo = alloc3(s, size, sizeof(HaloPeer3));
#ifdef ESF_MPI
    s->halo_requests = alloc3(s, 2 * size, sizeof(MPI_Request));
    s->halo_statuses = alloc3(s, 2 * size, sizeof(MPI_Status));
#endif
    int radius = 2 * c.grid.nx + 2, need_lo = s->first - radius,
        need_hi = s->first + s->count + radius;
    for (int peer = 0; peer < size; peer++) {
        int first = s->offsets[peer] / c.grid.ns, count = s->counts[peer] / c.grid.ns;
        s->coefficient_counts[peer] = 6 * count;
        s->coefficient_offsets[peer] = 6 * first;
        if (peer == rank)
            continue;
        int lo = first > need_lo ? first : need_lo,
            hi = first + count < need_hi ? first + count : need_hi;
        int other_lo = first - radius, other_hi = first + count + radius;
        int slo = s->first > other_lo ? s->first : other_lo,
            shi = s->first + s->count < other_hi ? s->first + s->count : other_hi;
        if (hi <= lo && shi <= slo)
            continue;
        HaloPeer3 *h = s->halo + s->halo_peers++;
        h->peer = peer;
        h->recv_first = lo;
        h->recv_count = hi > lo ? hi - lo : 0;
        h->send_first = slo;
        h->send_count = shi > slo ? shi - slo : 0;
        h->recv = alloc3(s, 2 * (size_t)h->recv_count * c.grid.ns, sizeof(double));
        h->send = alloc3(s, 2 * (size_t)h->send_count * c.grid.ns, sizeof(double));
    }
    s->rho = alloc3(s, s->total, sizeof(double));
    s->momentum = alloc3(s, s->total, sizeof(double));
    s->stage_n = alloc3(s, s->total, sizeof(double));
    s->stage_m = alloc3(s, s->total, sizeof(double));
    s->old_n = alloc3(s, s->total, sizeof(double));
    s->old_m = alloc3(s, s->total, sizeof(double));
    s->rhs_n = alloc3(s, s->total, sizeof(double));
    s->rhs_m = alloc3(s, s->total, sizeof(double));
    s->first_rhs_m = alloc3(s, s->total, sizeof(double));
    s->coeff = alloc3(s, 6 * columns, sizeof(double));
    s->local_coeff = alloc3(s, 6 * columns, sizeof(double));
    s->phi = alloc3(s, (c.grid.nx + 1) * (c.grid.ny + 1) + 1, sizeof(double));
    for (int a = 0; a < 3; a++) {
        int n = tube_face_count(&s->g, a);
        s->flux_n[a] = alloc3(s, n, sizeof(double));
        s->flux_m[a] = alloc3(s, n, sizeof(double));
        s->rate_face[a] = alloc3(s, n, sizeof(double));
    }
    if ((!rank || c.petsc_potential) &&
        (c.periodic_x ? potential_create_periodic(&s->p, c.grid.nx, c.grid.ny)
                      : potential_create(&s->p, c.grid.nx, c.grid.ny)))
        fail3(s, "potential allocation failed");
    /* Regional continuation also applies to gravity-driven polarization when
     * the imposed eastward field is zero (in particular an IGRF sector). */
    s->p.neumann_x = !c.periodic_x && c.regional_x;
#ifdef ESF_PETSC
    if (c.petsc_potential && solver3_petsc_create(s))
        fail3(s, "PETSc allocation failed");
#endif
    for (int p = s->first; p < s->first + s->count; p++) {
        int i = p % c.grid.nx, j = p / c.grid.nx;
        double apex = c.grid.hmin + (c.grid.hmax - c.grid.hmin) * (j + .5) / c.grid.ny;
        double longitudinal =
            cos(2 * PI * ((i + .5) / c.grid.nx - c.seed_offset_deg / c.grid.span));
        if (c.seed_east > 0) {
            double x = s->east_scale_m *
                       ((i + .5) / c.grid.nx - .5 - c.seed_offset_deg / c.grid.span) / c.seed_east;
            longitudinal = (2 * x * x - 1) * exp(-x * x);
        }
        double seed =
            c.amplitude * longitudinal * exp(-pow((apex - c.seed_height) / c.seed_width, 2));
        for (int k = 0; k < c.grid.ns; k++) {
            int id = p * c.grid.ns + k;
            s->rho[id] = solver3_initial_density(&c, s->g.cell[id].height) * (1 + seed);
        }
    }
    solver3_gather(s, s->rho);
    s->mass0 = solver3_mass(s, s->rho);
    if (!rank) {
        double lat, lon, h;
        geodetic_coordinates(s->target_apex, &lat, &lon, &h);
        printf(
            "Geometry ready: cells=%d, alignment=%.3g, E_parallel metric=%.3g; HK-350km apex h=%.2f km lon=%.3f\n",
            s->total, s->g.max_alignment_error, s->g.max_parallel_electric_metric_error,
            (vnorm(s->target_apex) - 6371200.) / 1000, lon);
        fflush(stdout);
    }
}
/**
 * Release uniform solver fields, geometry, halos, reconstruction weights and potential objects.
 */
void solver3_free(Solver3 *s) {
#ifdef ESF_PETSC
    if (s->petsc_context)
        solver3_petsc_destroy(s);
#endif
#ifdef ESF_MPI
    free(s->halo_requests);
    free(s->halo_statuses);
#endif
    tube_free(&s->g);
    potential_free(&s->p);
    free(s->counts);
    free(s->offsets);
    free(s->rho);
    free(s->momentum);
    free(s->stage_n);
    free(s->stage_m);
    free(s->old_n);
    free(s->old_m);
    free(s->rhs_n);
    free(s->rhs_m);
    free(s->first_rhs_m);
    free(s->coeff);
    free(s->local_coeff);
    free(s->phi);
    free(s->cached_velocity);
    free(s->cached_logrho);
    for (int k = 0; k < s->halo_peers; k++) {
        free(s->halo[k].recv);
        free(s->halo[k].send);
    }
    free(s->halo);
    free(s->coefficient_counts);
    free(s->coefficient_offsets);
    for (int a = 0; a < 3; a++) {
        free(s->flux_n[a]);
        free(s->flux_m[a]);
        free(s->rate_face[a]);
        free(s->pdm_weights[a]);
    }
}
/**
 * Collectively verify finite, nonnegative uniform density and consistent vacuum momentum.
 */
int solver3_valid(Solver3 *s, const double *n, const double *m) {
    double tick = wall3();
    int ok = 1;
    for (int k = s->first * s->c.grid.ns; k < (s->first + s->count) * s->c.grid.ns; k++)
        if (!isfinite(n[k]) || !isfinite(m[k]) || n[k] < 0 || (n[k] == 0 && m[k] != 0))
            ok = 0;
#ifdef ESF_MPI
    MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
    s->timing[PROF_VALID] += wall3() - tick;
    return ok;
}
/**
 * Volume-integrate uniform density and sum owned contributions across MPI ranks.
 */
double solver3_mass(Solver3 *s, const double *n) {
    double v = 0;
    for (int k = s->first * s->c.grid.ns; k < (s->first + s->count) * s->c.grid.ns; k++)
        v += n[k] * s->g.cell[k].volume;
#ifdef ESF_MPI
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    return v;
}
/**
 * Read a transverse potential node using this solver's boundary/periodic indexing convention.
 */
static double ph(const Solver3 *s, int i, int j) {
    return s->phi[j * (s->c.grid.nx + 1) + i];
}
/**
 * Compute transverse potential differences at a uniform geometric face with its boundary
 * convention.
 */
static void derivatives(const Solver3 *s, int a, int i, int j, double d[2]) {
    int nx = s->c.grid.nx, ny = s->c.grid.ny;
    if (a == 0) {
        int il = i ? i - 1 : i, ir = i < nx ? i + 1 : i, denom = ir - il;
        if (s->c.periodic_x && (i == 0 || i == nx)) {
            il = nx - 1;
            ir = 1;
            denom = 2;
            i = 0;
        }
        d[0] =
            nx * .5 * (ph(s, ir, j) + ph(s, ir, j + 1) - ph(s, il, j) - ph(s, il, j + 1)) / denom;
        d[1] = ny * (ph(s, i, j + 1) - ph(s, i, j));
    } else if (a == 1) {
        int jl = j ? j - 1 : j, jr = j < ny ? j + 1 : j;
        d[0] = nx * (ph(s, i + 1, j) - ph(s, i, j));
        d[1] = ny * .5 * (ph(s, i, jr) + ph(s, i + 1, jr) - ph(s, i, jl) - ph(s, i + 1, jl)) /
               (jr - jl);
    } else {
        d[0] = nx * .5 * (ph(s, i + 1, j) + ph(s, i + 1, j + 1) - ph(s, i, j) - ph(s, i, j + 1));
        d[1] = ny * .5 * (ph(s, i, j + 1) + ph(s, i + 1, j + 1) - ph(s, i, j) - ph(s, i + 1, j));
    }
}
/**
 * Return the sign-consistent slope with smallest magnitude; return zero at an extremum.
 */
static double minmod(double a, double b, double c) {
    return a * b > 0 && a * c > 0 ? copysign(fmin(fabs(a), fmin(fabs(b), fabs(c))), a) : 0;
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
 * Read a clamped or periodic scalar state sample.
 */
static double value(const Solver3 *s, const double *f, int i, int j, int k) {
    if (s->c.periodic_x) {
        if (i < 0 || i >= s->c.grid.nx)
            i = (i % s->c.grid.nx + s->c.grid.nx) % s->c.grid.nx;
    } else {
        if (i < 0)
            i = 0;
        if (i >= s->c.grid.nx)
            i = s->c.grid.nx - 1;
    }
    if (j < 0)
        j = 0;
    if (j >= s->c.grid.ny)
        j = s->c.grid.ny - 1;
    if (k < 0)
        k = 0;
    if (k >= s->c.grid.ns)
        k = s->c.grid.ns - 1;
    return f[tube_cell_id(&s->g, i, j, k)];
}
/**
 * Reconstruct a uniform face state using the selected slope policy.
 */
static double state(const Solver3 *s, const double *f, int a, int i, int j, int k, int side) {
    int p[3] = {i, j, k}, l[3] = {i, j, k}, r[3] = {i, j, k};
    l[a]--;
    r[a]++;
    double x = value(s, f, p[0], p[1], p[2]), u = value(s, f, l[0], l[1], l[2]),
           v = value(s, f, r[0], r[1], r[2]);
    return x + .5 * side * minmod(1.3 * (x - u), .5 * (v - u), 1.3 * (v - x));
}
/**
 * Reconstruct a bounded transverse density face value with physical-volume PDM4 or the
 * configured lower-order stencil.
 */
static double exb_state(const Solver3 *s, const double *f, int axis, int i, int j, int k, int fid,
                        int side) {
    int donor[3] = {i, j, k};
    donor[axis] -= side > 0;
    if (!s->c.exb_pdm4 || !isfinite(s->pdm_weights[axis][4 * fid]))
        return state(s, f, axis, donor[0], donor[1], donor[2], side);
    const double *w = s->pdm_weights[axis] + 4 * fid;
    double high = 0;
    for (int c = 0; c < 4; c++) {
        int cell[3] = {i, j, k};
        cell[axis] += c - 2;
        high += w[c] * value(s, f, cell[0], cell[1], cell[2]);
    }
    int before[3] = {donor[0], donor[1], donor[2]}, after[3] = {donor[0], donor[1], donor[2]};
    before[axis] -= side;
    after[axis] += side;
    return tube_pdm_limit(high, value(s, f, before[0], before[1], before[2]),
                          value(s, f, donor[0], donor[1], donor[2]),
                          value(s, f, after[0], after[1], after[2]), s->c.pdm_beta);
}
/**
 * Return parallel velocity with an exact-vacuum zero branch.
 */
static double velocity_value(const Solver3 *s, const double *n, const double *m, int i, int j,
                             int k) {
    (void)n;
    (void)m;
    return value(s, s->cached_velocity, i, j, k);
}
/**
 * Construct the bounded parallel-velocity face state used for momentum transport.
 */
static double velocity_state(const Solver3 *s, const double *n, const double *m, int a, int i,
                             int j, int k, int side) {
    int l[3] = {i, j, k}, r[3] = {i, j, k};
    l[a]--;
    r[a]++;
    double x = velocity_value(s, n, m, i, j, k), u = velocity_value(s, n, m, l[0], l[1], l[2]),
           v = velocity_value(s, n, m, r[0], r[1], r[2]);
    return x + .5 * side * minmod(1.3 * (x - u), .5 * (v - u), 1.3 * (v - x));
}
/**
 * Read the physical length of an along-field segment.
 */
static double line_length(const Solver3 *s, int i, int j, int k) {
    double a = s->g.face[2][tube_face_id(&s->g, 2, i, j, k)].parallel_area,
           b = s->g.face[2][tube_face_id(&s->g, 2, i, j, k + 1)].parallel_area;
    return s->g.cell[tube_cell_id(&s->g, i, j, k)].volume / (.5 * (fabs(a) + fabs(b)));
}
/**
 * Compute a positive-state logarithmic pressure jump or its exact-vacuum secant extension.
 */
static double log_difference(double a, double b, double local, double la, double lb) {
    return a > 0 && b > 0 ? lb - la : (b - a) / local;
}
/**
 * Construct a nonuniform parallel logarithmic density gradient for the uniform isothermal
 * pressure force.
 */
static double pressure_log_gradient(const Solver3 *s, const double *n, int i, int j, int k) {
    int ns = s->c.grid.ns, lo = k ? k - 1 : k, hi = k + 1 < ns ? k + 1 : k;
    double x = n[tube_cell_id(&s->g, i, j, k)];
    if (x == 0)
        return 0;
    double own = line_length(s, i, j, k), hl = .5 * (own + line_length(s, i, j, lo)),
           hr = .5 * (own + line_length(s, i, j, hi));
    int il = tube_cell_id(&s->g, i, j, lo), ic = tube_cell_id(&s->g, i, j, k),
        ir = tube_cell_id(&s->g, i, j, hi);
    double dl = log_difference(n[il], x, x, s->cached_logrho[il], s->cached_logrho[ic]),
           dr = log_difference(x, n[ir], x, s->cached_logrho[ic], s->cached_logrho[ir]);
    if (k == 0)
        return dr / hr;
    if (k == ns - 1)
        return dl / hl;
    return (hr / hl * dl + hl / hr * dr) / (hl + hr);
}
/**
 * Exchange current state, integrate conductivity, solve the potential, and evaluate the uniform
 * transport/source RHS.
 */
int solver3_rhs(Solver3 *s, double *n, double *m, double *dn, double *dm) {
    double tick = wall3();
    exchange_state(s, n, m);
    s->timing[PROF_STATE] += wall3() - tick;
    tick = wall3();
    potential_coefficients(&s->g, n, s->first, s->count, s->c.hall, s->local_coeff);
    s->timing[PROF_COEFF] += wall3() - tick;
    tick = wall3();
#ifdef ESF_MPI
    if (s->c.petsc_potential)
        MPI_Allgatherv(s->local_coeff + 6 * s->first, 6 * s->count, MPI_DOUBLE, s->coeff,
                       s->coefficient_counts, s->coefficient_offsets, MPI_DOUBLE, MPI_COMM_WORLD);
    else
        MPI_Gatherv(s->local_coeff + 6 * s->first, 6 * s->count, MPI_DOUBLE, s->coeff,
                    s->coefficient_counts, s->coefficient_offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    memcpy(s->coeff, s->local_coeff, 6 * s->c.grid.nx * s->c.grid.ny * sizeof(double));
#endif
    s->timing[PROF_COEFF_MPI] += wall3() - tick;
#ifdef ESF_PETSC
    if (s->c.petsc_potential) {
        if (solver3_petsc_solve(s))
            return -1;
        s->solves++;
        return solver3_transport(s, n, m, dn, dm);
    }
#endif
    int code = 0;
    if (!s->rank) {
        tick = wall3();
        /* phi_total=phi_polarization+phi_background. Move the full tensor
         * background current to the RHS before solving polarization. */
        for (int p = 0; p < s->c.grid.nx * s->c.grid.ny; p++) {
            s->coeff[6 * p + 4] -= s->coeff[6 * p] * s->background_potential_gradient;
            s->coeff[6 * p + 5] -= s->coeff[6 * p + 2] * s->background_potential_gradient;
        }
        potential_assemble(&s->p, s->coeff);
        s->timing[PROF_ASSEMBLY] += wall3() - tick;
        tick = wall3();
        code = potential_solve(&s->p, s->c.rtol);
        s->timing[PROF_SOLVE] += wall3() - tick;
        s->potential_iterations += s->p.iterations;
        if (code >= 0)
            for (int p = 0; p < s->p.n; p++)
                s->phi[p] =
                    s->p.phi[p] + s->background_potential_gradient *
                                      ((double)(p % (s->c.grid.nx + 1)) / s->c.grid.nx - .5);
        s->max_potential_residual = fmax(s->max_potential_residual, s->p.residual);
    }
    tick = wall3();
    int phi_n = (s->c.grid.nx + 1) * (s->c.grid.ny + 1);
    if (!s->rank)
        s->phi[phi_n] = code;
    broadcast(s, s->phi, (phi_n + 1) * sizeof(double));
    code = (int)s->phi[phi_n];
    if (code < 0)
        return -1;
    s->solves++;
    s->timing[PROF_BCAST] += wall3() - tick;
    return solver3_transport(s, n, m, dn, dm);
}
/**
 * Evaluate conservative uniform particle/momentum face fluxes, gravity/pressure sources,
 * boundary budgets and CFL rates.
 */
int solver3_transport(Solver3 *s, const double *n, const double *m, double *dn, double *dm) {
    double tick = wall3();
    double cs = sqrt(1.380649e-23 * s->c.temperature / (16 * 1.66053906660e-27));
    int nx = s->c.grid.nx, ny = s->c.grid.ny, ns = s->c.grid.ns;
    if (s->c.parallel_relaxation) {
        if (!s->cached_velocity)
            s->cached_velocity = alloc3(s, s->total, sizeof(double));
        if (!s->cached_logrho)
            s->cached_logrho = alloc3(s, s->total, sizeof(double));
        /* Every owned adjacent face reconstructs at most two columns away.
         * Include complete +/-2-row halos, also covering periodic x wraps. */
        int begin = s->first - 2 * nx - 2, end = s->first + s->count + 2 * nx + 2;
        if (begin < 0)
            begin = 0;
        if (end > nx * ny)
            end = nx * ny;
        for (int k = begin * ns; k < end * ns; k++) {
            s->cached_velocity[k] = n[k] > 0 ? m[k] / n[k] : 0;
            s->cached_logrho[k] = n[k] > 0 ? log(n[k]) : 0;
        }
    }
    for (int a = 0; a < 3; a++)
        for (int j = 0; j < ny + (a == 1); j++)
            for (int i = 0; i < nx + (a == 0); i++) {
                int left = -1, right = -1;
                if (a == 0) {
                    if (i)
                        left = j * nx + i - 1;
                    if (i < nx)
                        right = j * nx + i;
                    if (s->c.periodic_x) {
                        if (!i)
                            left = j * nx + nx - 1;
                        if (i == nx)
                            right = j * nx;
                    }
                } else if (a == 1) {
                    if (j)
                        left = (j - 1) * nx + i;
                    if (j < ny)
                        right = j * nx + i;
                } else
                    left = right = j * nx + i;
                if (!(left >= s->first && left < s->first + s->count) &&
                    !(right >= s->first && right < s->first + s->count))
                    continue;
                double d[2];
                derivatives(s, a, i, j, d);
                for (int k = 0; k < ns + (a == 2); k++) {
                    int p[3] = {i, j, k}, l[3] = {i, j, k};
                    l[a]--;
                    int fid = tube_face_id(&s->g, a, i, j, k);
                    /* Both periodic copies use the same geometric face coefficients,
                     * so their fluxes cancel exactly, including across MPI owners. */
                    int geometry_id = a == 0 && s->c.periodic_x && i == nx
                                          ? tube_face_id(&s->g, a, 0, j, k)
                                          : fid;
                    const TubeFace *f = s->g.face[a] + geometry_id;
                    double q = f->e[0] * d[0] + f->e[1] * d[1];
                    double nl = value(s, n, l[0], l[1], l[2]), nr = value(s, n, p[0], p[1], p[2]);
                    double ml = value(s, m, l[0], l[1], l[2]), mr = value(s, m, p[0], p[1], p[2]);
                    double ul = nl > 0 ? ml / nl : 0, ur = nr > 0 ? mr / nr : 0;
                    int wall = (a == 0 && !s->c.periodic_x && (i == 0 || i == nx)) ||
                               (a == 1 && (j == 0 || j == ny)) || (a == 2 && (k == 0 || k == ns));
                    int open = (a == 1 && (s->c.background_east != 0 || s->c.open_apex)) ||
                               (a == 0 && !s->c.periodic_x && s->c.regional_x) ||
                               (a == 2 && s->c.open_ends);
                    if (wall && !open)
                        q = 0;
                    double flow_n = 0, flow_m = 0;
                    if (q != 0) {
                        int side = q > 0 ? 1 : -1;
                        flow_n = q * exb_state(s, n, a, i, j, k, geometry_id, side);
                        flow_m = s->c.parallel_relaxation
                                     ? flow_n * exb_state(s, s->cached_velocity, a, i, j, k,
                                                          geometry_id, side)
                                     : q * exb_state(s, m, a, i, j, k, geometry_id, side);
                    }
                    double speed = 0, parallel_rate = 0;
                    if (a == 2 && s->c.parallel) {
                        if (s->c.parallel_order == 2 && k > 0 && k < ns) {
                            nl = state(s, n, a, l[0], l[1], l[2], 1);
                            nr = state(s, n, a, p[0], p[1], p[2], -1);
                            ml = state(s, m, a, l[0], l[1], l[2], 1);
                            mr = state(s, m, a, p[0], p[1], p[2], -1);
                            ul = nl > 0 ? ml / nl : 0;
                            ur = nr > 0 ? mr / nr : 0;
                        }
                        if (s->c.parallel_relaxation && s->c.parallel_order == 2 && k > 0 &&
                            k < ns) {
                            ul = velocity_state(s, n, m, a, l[0], l[1], l[2], 1);
                            ur = velocity_state(s, n, m, a, p[0], p[1], p[2], -1);
                            ml = nl * ul;
                            mr = nr * ur;
                        }
                        if (k == 0 && !open) {
                            nl = nr;
                            ml = -mr;
                            ul = -ur;
                        }
                        if (k == ns && !open) {
                            nr = nl;
                            mr = -ml;
                            ur = -ul;
                        }
                        speed = fmax(fabs(ul), fabs(ur)) + cs;
                        double area = f->parallel_area;
                        double parallel_n =
                            .5 * area * (ml + mr) - .5 * fabs(area) * speed * (nr - nl);
                        if (!s->c.parallel_relaxation)
                            flow_m += .5 * area * (ml * ul + mr * ur + cs * cs * (nl + nr)) -
                                      .5 * fabs(area) * speed * (mr - ml);
                        parallel_rate = fabs(area) * speed;
                        if (s->c.parallel_relaxation && (!wall || open)) {
                            int kl = k ? k - 1 : 0, kr = k < ns ? k : ns - 1;
                            int il = tube_cell_id(&s->g, i, j, kl),
                                ir = tube_cell_id(&s->g, i, j, kr);
                            const TubeCell *cl = s->g.cell + il, *cr = s->g.cell + ir;
                            double nu = .5 * (cl->nu + cr->nu);
                            if (nu > 0) {
                                double al =
                                    .5 * (fabs(s->g.face[2][tube_face_id(&s->g, 2, i, j, kl)]
                                                   .parallel_area) +
                                          fabs(s->g.face[2][tube_face_id(&s->g, 2, i, j, kl + 1)]
                                                   .parallel_area));
                                double ar =
                                    .5 * (fabs(s->g.face[2][tube_face_id(&s->g, 2, i, j, kr)]
                                                   .parallel_area) +
                                          fabs(s->g.face[2][tube_face_id(&s->g, 2, i, j, kr + 1)]
                                                   .parallel_area));
                                double ds = .5 * (cl->volume / al + cr->volume / ar),
                                       g = .5 * (cl->gravity_parallel + cr->gravity_parallel);
                                double kn = nu * ds / cs, theta = 1 / (1 + kn * kn),
                                       D = cs * cs / nu, pe = g * ds / (cs * cs);
                                double bl = bernoulli(-pe), br = bernoulli(pe),
                                       conductance = fabs(area) * D / ds;
                                /* Strong drag: conservative Scharfetter-Gummel flux
                                 * rho*u = (rho*g - cs^2*d_s rho)/nu.
                                 * The quadratic blend changes the resolved Euler flux
                                 * by O(ds^2), and eliminates acoustic mass diffusion
                                 * as nu*ds/cs -> infinity. Both faces share this flux. */
                                double terminal = conductance * (bl * n[il] - br * n[ir]);
                                parallel_n = theta * parallel_n + (1 - theta) * terminal;
                                parallel_rate = fmax(parallel_rate,
                                                     theta * parallel_rate +
                                                         (1 - theta) * conductance * fmax(bl, br));
                            }
                        }
                        flow_n += parallel_n;
                        if (s->c.parallel_relaxation && (!wall || open))
                            flow_m += parallel_n * (parallel_n >= 0 ? ul : ur);
                        if (wall && !open)
                            flow_n = 0;
                    }
                    s->flux_n[a][fid] = flow_n;
                    s->flux_m[a][fid] = flow_m;
                    s->rate_face[a][fid] = fabs(q) + parallel_rate;
                }
            }
    s->timing[PROF_FACES] += wall3() - tick;
    tick = wall3();
    double rate = 0;
    for (int p = s->first; p < s->first + s->count; p++) {
        int i = p % nx, j = p / nx;
        for (int k = 0; k < ns; k++) {
            int id = p * ns + k;
            double rn = 0, rm = 0, local = 0;
            const TubeCell *c = s->g.cell + id;
            for (int a = 0; a < 3; a++) {
                int lo = tube_face_id(&s->g, a, i, j, k),
                    hi = tube_face_id(&s->g, a, i + (a == 0), j + (a == 1), k + (a == 2));
                rn += s->flux_n[a][lo] - s->flux_n[a][hi];
                rm += s->flux_m[a][lo] - s->flux_m[a][hi];
                local += s->rate_face[a][lo] + s->rate_face[a][hi];
            }
            dn[id] = rn / c->volume;
            dm[id] = rm / c->volume;
            if (s->c.parallel) {
                int z0 = tube_face_id(&s->g, 2, i, j, k), z1 = tube_face_id(&s->g, 2, i, j, k + 1);
                if (s->c.parallel_relaxation)
                    dm[id] += n[id] * (c->gravity_parallel -
                                       cs * cs * pressure_log_gradient(s, n, i, j, k));
                else
                    dm[id] +=
                        cs * cs * n[id] *
                            (s->g.face[2][z1].parallel_area - s->g.face[2][z0].parallel_area) /
                            c->volume +
                        n[id] * c->gravity_parallel;
            } else
                dm[id] = 0;
            rate = fmax(rate, local / c->volume);
        }
    }
    s->timing[PROF_CELLS] += wall3() - tick;
    tick = wall3();
#ifdef ESF_MPI
    MPI_Allreduce(MPI_IN_PLACE, &rate, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    s->boundary_rate = 0;
    if (s->c.background_east != 0 || s->c.open_apex || s->c.open_ends ||
        (!s->c.periodic_x && s->c.regional_x)) {
        if (s->c.background_east != 0 || s->c.open_apex) {
            for (int i = 0; i < nx; i++)
                for (int k = 0; k < ns; k++) {
                    if (i >= s->first && i < s->first + s->count)
                        s->boundary_rate += s->flux_n[1][tube_face_id(&s->g, 1, i, 0, k)];
                    int top = (ny - 1) * nx + i;
                    if (top >= s->first && top < s->first + s->count)
                        s->boundary_rate -= s->flux_n[1][tube_face_id(&s->g, 1, i, ny, k)];
                }
        }
        if (s->c.regional_x && !s->c.periodic_x)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < ns; k++) {
                    int left = j * nx, right = j * nx + nx - 1;
                    if (left >= s->first && left < s->first + s->count)
                        s->boundary_rate += s->flux_n[0][tube_face_id(&s->g, 0, 0, j, k)];
                    if (right >= s->first && right < s->first + s->count)
                        s->boundary_rate -= s->flux_n[0][tube_face_id(&s->g, 0, nx, j, k)];
                }
        if (s->c.open_ends)
            for (int p = s->first; p < s->first + s->count; p++) {
                int i = p % nx, j = p / nx;
                s->boundary_rate += s->flux_n[2][tube_face_id(&s->g, 2, i, j, 0)] -
                                    s->flux_n[2][tube_face_id(&s->g, 2, i, j, ns)];
            }
#ifdef ESF_MPI
        MPI_Allreduce(MPI_IN_PLACE, &s->boundary_rate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    }
    s->cfl_rate = rate;
    s->timing[PROF_REDUCE] += wall3() - tick;
    return 0;
}
