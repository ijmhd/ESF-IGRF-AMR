/*
 * Traced IGRF flux-tube hexahedra in ECEF metres. Apex longitude and apex spherical altitude
 * label a tube; the along-field coordinate is a weighted tracing distance, not geographic
 * latitude.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "tube.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define PI 3.14159265358979323846
#define R0 6371200.
/**
 * Return the unit magnetic-field tangent in ECEF for field-line integration.
 */
static V3 tangent(const Igrf *m, V3 x) {
    return vunit(igrf_field(m, x));
}
/**
 * Advance one signed physical-distance field-line step with classical RK4 using unit field
 * tangents.
 */
static V3 rk4(const Igrf *m, V3 x, double ds) {
    V3 a = tangent(m, x), b = tangent(m, vadd(x, vmul(a, .5 * ds)));
    V3 c = tangent(m, vadd(x, vmul(b, .5 * ds))), d = tangent(m, vadd(x, vmul(c, ds)));
    return vadd(x, vmul(vadd(vadd(a, vmul(b, 2)), vadd(vmul(c, 2), d)), ds / 6));
}
/**
 * Return B dot radial-unit-vector; its zero locates a spherical-altitude apex.
 */
static double radial(const Igrf *m, V3 x) {
    return vdot(igrf_field(m, x), vunit(x));
}
/**
 * Trace from a given position toward increasing radius and bracket/bisect the radial-field zero;
 * reject escaping or interior trajectories.
 */
int tube_apex(const Igrf *m, V3 start, V3 *apex) {
    double direction = radial(m, start) < 0 ? -1 : 1;
    V3 x = start;
    for (int i = 0; i < 10000; i++) {
        V3 y = rk4(m, x, 5000 * direction);
        if (radial(m, x) * radial(m, y) <= 0) {
            double lo = 0, hi = 1;
            for (int j = 0; j < 45; j++) {
                double t = .5 * (lo + hi);
                V3 p = vadd(x, vmul(vsub(y, x), t));
                if (radial(m, p) * radial(m, x) > 0)
                    lo = t;
                else
                    hi = t;
            }
            *apex = vadd(x, vmul(vsub(y, x), .5 * (lo + hi)));
            return 0;
        }
        x = y;
        if (vnorm(x) > R0 + 1e7 || vnorm(x) < R0)
            return -1;
    }
    return -1;
}
typedef struct {
    V3 x;
    double w;
} TracePoint;
/**
 * Trace one conjugate half tube from its apex to the spherical base boundary, accumulating the
 * altitude-weighted along-field coordinate.
 */
static int half_trace(const Igrf *m, V3 x, int sign, double base, double step, TracePoint *t) {
    int n = 1;
    t[0] = (TracePoint){x, 0};
    while (vnorm(x) > R0 + base && n < 4096) {
        double h = vnorm(x) - R0, ds = fmin(step, fmax(1000., .05 * (h - base) + 1000.));
        V3 y = rk4(m, x, sign * ds);
        if (vnorm(y) < R0 + base) {
            double lo = 0, hi = 1;
            for (int j = 0; j < 40; j++) {
                double a = .5 * (lo + hi);
                if (vnorm(vadd(x, vmul(vsub(y, x), a))) > R0 + base)
                    lo = a;
                else
                    hi = a;
            }
            y = vadd(x, vmul(vsub(y, x), .5 * (lo + hi)));
        }
        double hh = .5 * (vnorm(x) + vnorm(y)) - R0;
        double scale = fmax(8000., fmin(80000., 8000. + .12 * (hh - base)));
        t[n] = (TracePoint){y, t[n - 1].w + vnorm(vsub(y, x)) / scale};
        x = y;
        n++;
        if (vnorm(x) <= R0 + base + 1e-5)
            break;
    }
    return n == 4096 ? -1 : n;
}
/**
 * Find the magnetic-apex latitude at specified geographic apex longitude and spherical altitude
 * by a radial-field root solve.
 */
int tube_apex_label(const Igrf *m, double lon, double height, V3 *apex) {
    double lo = -45, hi = 45, r = R0 + height;
    if (radial(m, spherical_position(lo, lon, r)) * radial(m, spherical_position(hi, lon, r)) >= 0)
        return -1;
    for (int n = 0; n < 50; n++) {
        double a = .5 * (lo + hi);
        if (radial(m, spherical_position(a, lon, r)) > 0)
            lo = a;
        else
            hi = a;
    }
    *apex = spherical_position(.5 * (lo + hi), lon, r);
    return 0;
}
/**
 * Trace both conjugate halves and resample ns+1 nodes uniformly in the weighted coordinate;
 * output spans the two base-altitude endpoints.
 */
int tube_trace(const Igrf *m, double lon, double height, double base, double step, int ns,
               V3 *out) {
    V3 apex;
    if (tube_apex_label(m, lon, height, &apex))
        return -1;
    TracePoint *south = calloc(4096, sizeof(*south)), *north = calloc(4096, sizeof(*north));
    if (!south || !north) {
        free(south);
        free(north);
        return -2;
    }
    int sn = half_trace(m, apex, -1, base, step, south),
        nn = half_trace(m, apex, 1, base, step, north);
    if (sn < 0 || nn < 0) {
        free(south);
        free(north);
        return -3;
    }
    double sw = south[sn - 1].w, nw = north[nn - 1].w;
    for (int k = 0; k <= ns; k++) {
        double w = (sw + nw) * k / ns;
        TracePoint *a;
        int n;
        if (w < sw) {
            w = sw - w;
            a = south;
            n = sn;
        } else {
            w -= sw;
            a = north;
            n = nn;
        }
        int j = 0;
        while (j + 1 < n && a[j + 1].w < w)
            j++;
        if (j + 1 == n)
            out[k] = a[j].x;
        else {
            double f = (w - a[j].w) / (a[j + 1].w - a[j].w);
            out[k] = vadd(a[j].x, vmul(vsub(a[j + 1].x, a[j].x), f));
        }
    }
    free(south);
    free(north);
    return 0;
}
/**
 * Flatten node coordinates with the along-field coordinate contiguous.
 */
static int nid(const TubeGrid *g, int i, int j, int k) {
    return ((j * (g->c.nx + 1) + i) * (g->c.ns + 1) + k);
}
/**
 * Flatten cell coordinates in [ny,nx,ns] storage with k contiguous.
 */
int tube_cell_id(const TubeGrid *g, int i, int j, int k) {
    return ((j * g->c.nx + i) * g->c.ns + k);
}
/**
 * Flatten an axis-specific staggered face index, accounting for the extra face in its normal
 * direction.
 */
int tube_face_id(const TubeGrid *g, int axis, int i, int j, int k) {
    return (j * (g->c.nx + (axis == 0)) + i) * (g->c.ns + (axis == 2)) + k;
}
/**
 * Return the number of faces in one staggered axis family.
 */
int tube_face_count(const TubeGrid *g, int a) {
    return (g->c.nx + (a == 0)) * (g->c.ny + (a == 1)) * (g->c.ns + (a == 2));
}
/**
 * Differentiate the trilinear eight-node cell mapping to obtain its three covariant basis
 * vectors.
 */
static void jacobian(const V3 v[8], double x, double y, double z, V3 e[3]) {
    e[0] = e[1] = e[2] = (V3){0, 0, 0};
    for (int a = 0; a < 2; a++)
        for (int b = 0; b < 2; b++)
            for (int c = 0; c < 2; c++) {
                V3 p = v[(a * 2 + b) * 2 + c];
                double wx = a ? x : 1 - x, wy = b ? y : 1 - y, wz = c ? z : 1 - z;
                e[0] = vadd(e[0], vmul(p, (a ? 1 : -1) * wy * wz));
                e[1] = vadd(e[1], vmul(p, wx * (b ? 1 : -1) * wz));
                e[2] = vadd(e[2], vmul(p, wx * wy * (c ? 1 : -1)));
            }
}
/**
 * Average the four corners of a staggered face in physical ECEF coordinates.
 */
static V3 face_center(const TubeGrid *g, int a, int p[3]) {
    int b = (a + 1) % 3, c = (a + 2) % 3;
    V3 x = {0, 0, 0};
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            int q[3] = {p[0], p[1], p[2]};
            q[b] += i;
            q[c] += j;
            x = vadd(x, g->node[nid(g, q[0], q[1], q[2])]);
        }
    return vmul(x, .25);
}
/**
 * Allocate nodes/cells/faces for a grid and initialize metric diagnostics; no field tracing is
 * performed here.
 */
int tube_geometry_allocate(TubeGrid *g, GridConfig c) {
    memset(g, 0, sizeof(*g));
    g->c = c;
    g->min_volume = HUGE_VAL;
    g->node = calloc((size_t)(c.nx + 1) * (c.ny + 1) * (c.ns + 1), sizeof(V3));
    g->cell = calloc((size_t)c.nx * c.ny * c.ns, sizeof(TubeCell));
    for (int a = 0; a < 3; a++)
        g->face[a] = calloc(tube_face_count(g, a), sizeof(TubeFace));
    if (!g->node || !g->cell || !g->face[0] || !g->face[1] || !g->face[2])
        return -1;
    return 0;
}
/**
 * Trace every node line then construct cell/face metrics over the full grid.
 */
int tube_geometry(TubeGrid *g, const Igrf *m, GridConfig c) {
    if (tube_geometry_allocate(g, c))
        return -1;
    for (int j = 0; j <= c.ny; j++)
        for (int i = 0; i <= c.nx; i++)
            if (tube_trace(m, c.lon + c.span * ((double)i / c.nx - .5),
                           c.hmin + (c.hmax - c.hmin) * j / c.ny, c.hbase, c.trace_step, c.ns,
                           g->node + nid(g, i, j, 0)))
                return -2;
    return tube_geometry_rows(g, m, 0, c.ny);
}
/**
 * Build metrics for a disjoint cell-row range: Jacobian-volume quadrature, transverse covectors,
 * field-aligned tangent, gravity and integrated face flux factors. Reject folded/degenerate
 * geometry.
 */
int tube_geometry_rows(TubeGrid *g, const Igrf *m, int first, int last) {
    GridConfig c = g->c;
    if (first < 0 || last > c.ny || first >= last)
        return -7;
    const double q[2] = {.5 - .5 / sqrt(3.), .5 + .5 / sqrt(3.)};
    for (int j = first; j < last; j++)
        for (int i = 0; i < c.nx; i++)
            for (int k = 0; k < c.ns; k++) {
                V3 v[8], x = {0, 0, 0}, e[3];
                for (int a = 0; a < 2; a++)
                    for (int b = 0; b < 2; b++)
                        for (int d = 0; d < 2; d++) {
                            v[(a * 2 + b) * 2 + d] = g->node[nid(g, i + a, j + b, k + d)];
                            x = vadd(x, vmul(v[(a * 2 + b) * 2 + d], .125));
                        }
                TubeCell *cell = g->cell + tube_cell_id(g, i, j, k);
                cell->position = x;
                jacobian(v, .5, .5, .5, e);
                double det = vdot(e[0], vcross(e[1], e[2]));
                if (!isfinite(det) || det == 0)
                    return -3;
                cell->gx = vmul(vcross(e[1], e[2]), 1. / (det * c.nx));
                cell->gy = vmul(vcross(e[2], e[0]), 1. / (det * c.ny));
                V3 field = igrf_field(m, x), bf = vunit(field), b = vunit(e[2]);
                if (vdot(b, bf) < 0)
                    b = vmul(b, -1);
                cell->b = b;
                cell->B = vnorm(field);
                cell->height = vnorm(x) - R0;
                cell->gravity_parallel = vdot(vmul(vunit(x), -9.8 * pow(R0 / vnorm(x), 2)), b);
                g->max_alignment_error = fmax(g->max_alignment_error, vnorm(vcross(b, bf)));
                g->max_parallel_electric_metric_error =
                    fmax(g->max_parallel_electric_metric_error,
                         fmax(fabs(vdot(b, cell->gx)) / vnorm(cell->gx),
                              fabs(vdot(b, cell->gy)) / vnorm(cell->gy)));
                for (int a = 0; a < 2; a++)
                    for (int b0 = 0; b0 < 2; b0++)
                        for (int d = 0; d < 2; d++) {
                            jacobian(v, q[a], q[b0], q[d], e);
                            double determinant = vdot(e[0], vcross(e[1], e[2]));
                            if (det * determinant <= 0)
                                return -4;
                            cell->volume += fabs(determinant) / 8.;
                        }
                g->min_volume = fmin(g->min_volume, cell->volume);
            }
    int n[3] = {c.nx, c.ny, c.ns};
    for (int a = 0; a < 3; a++)
        for (int j = first; j < last + (a == 1 && last == c.ny); j++)
            for (int i = 0; i < c.nx + (a == 0); i++)
                for (int k = 0; k < c.ns + (a == 2); k++) {
                    int p[3] = {i, j, k}, b = (a + 1) % 3, d = (a + 2) % 3;
                    V3 e[3] = {{0, 0, 0}}, x = face_center(g, a, p);
                    for (int t = 0; t < 2; t++) {
                        int left[3] = {i, j, k}, right[3] = {i, j, k};
                        left[d] += t;
                        right[d] += t;
                        right[b]++;
                        e[b] = vadd(e[b], vmul(vsub(g->node[nid(g, right[0], right[1], right[2])],
                                                    g->node[nid(g, left[0], left[1], left[2])]),
                                               .5 * n[b]));
                        left[0] = right[0] = i;
                        left[1] = right[1] = j;
                        left[2] = right[2] = k;
                        left[b] += t;
                        right[b] += t;
                        right[d]++;
                        e[d] = vadd(e[d], vmul(vsub(g->node[nid(g, right[0], right[1], right[2])],
                                                    g->node[nid(g, left[0], left[1], left[2])]),
                                               .5 * n[d]));
                    }
                    int pp[3] = {i, j, k}, pm[3] = {i, j, k};
                    pp[a] = p[a] < n[a] ? p[a] + 1 : p[a];
                    pm[a] = p[a] > 0 ? p[a] - 1 : p[a];
                    e[a] = vmul(vsub(face_center(g, a, pp), face_center(g, a, pm)),
                                (double)n[a] / (pp[a] - pm[a]));
                    double det = vdot(e[0], vcross(e[1], e[2]));
                    if (det == 0 || !isfinite(det))
                        return -5;
                    V3 gx = vmul(vcross(e[1], e[2]), 1 / det),
                       gy = vmul(vcross(e[2], e[0]), 1 / det);
                    V3 normal = vmul(vcross(e[b], e[d]), (det > 0 ? 1. : -1.) / (n[b] * n[d]));
                    V3 bf = igrf_field(m, x), field = vunit(e[2]);
                    if (vdot(field, bf) < 0)
                        field = vmul(field, -1);
                    TubeFace *f = g->face[a] + tube_face_id(g, a, i, j, k);
                    f->e[0] = -vdot(normal, vcross(gx, field)) / vnorm(bf);
                    f->e[1] = -vdot(normal, vcross(gy, field)) / vnorm(bf);
                    f->parallel_area = a == 2 ? vdot(normal, field) : 0.;
                    f->magnetic_flux = f->parallel_area * vnorm(bf);
                }
    for (int j = first; j < last; j++)
        for (int i = 0; i < c.nx; i++)
            for (int k = 0; k < c.ns; k++) {
                double a = g->face[2][tube_face_id(g, 2, i, j, k)].magnetic_flux,
                       b = g->face[2][tube_face_id(g, 2, i, j, k + 1)].magnetic_flux;
                if (a <= 0 || b <= 0)
                    return -6;
                g->max_magnetic_flux_jump =
                    fmax(g->max_magnetic_flux_jump, fabs(b - a) / fmax(a, b));
            }
    return 0;
}
/**
 * Free every canonical geometric array and clear the grid object.
 */
void tube_free(TubeGrid *g) {
    free(g->node);
    free(g->cell);
    for (int a = 0; a < 3; a++)
        free(g->face[a]);
    memset(g, 0, sizeof(*g));
}
/**
 * Evaluate the code's normalized Chapman-like layer exp(1-z-exp(-z)), with z=(height-
 * peak)/scale.
 */
double chapman(double height, double peak, double scale) {
    double z = (height - peak) / scale;
    return exp(1 - z - exp(-z));
}
/**
 * Solve the steady perpendicular collisional gravity drift for a species with signed
 * gyrofrequency.
 */
static V3 gravity_drift(V3 g, V3 b, double nu, double omega) {
    g = vsub(g, vmul(b, vdot(g, b)));
    return vmul(vsub(vmul(g, nu), vmul(vcross(b, g), omega)), 1 / (nu * nu + omega * omega));
}
/**
 * Build O+/fixed-molecular Pedersen/Hall and gravity-current responses, including electrons
 * under quasineutrality.
 */
static void conductivity_cell(TubeCell *c, double e_peak, double nuF, double nue) {
    const double qe = 1.602176634e-19, amu = 1.66053906660e-27, me = 9.1093837139e-31;
    double h = c->height;
    double nuE = 600 * exp(-(h - 110e3) / 10e3) + 1.2 * 26.7 * exp(-h / 81400.);
    double oi = qe * c->B / (16 * amu), om = qe * c->B / (30 * amu), oe = qe * c->B / me;
    double ki = oi / nuF, km = om / nuE, ke = oe / nue;
    double p_i = ki / (1 + ki * ki), p_m = km / (1 + km * km), p_e = ke / (1 + ke * ke);
    /* Convention J_H = sigma_H b cross E. Quasineutral form avoids
     * cancellation between magnetized ion/electron Hall responses. */
    double hi = 1 / (1 + ki * ki) - 1 / (1 + ke * ke), hm = 1 / (1 + km * km) - 1 / (1 + ke * ke);
    c->nu = nuF;
    c->nE = e_peak * chapman(h, 110e3, 12e3);
    c->pF = qe / c->B * 1e12 * (p_i + p_e);
    c->hF = qe / c->B * 1e12 * hi;
    c->pE = qe / c->B * c->nE * (p_m + p_e);
    c->hE = qe / c->B * c->nE * hm;
    V3 grav = vmul(vunit(c->position), -9.8 * pow(R0 / vnorm(c->position), 2));
    V3 jF = vmul(vsub(gravity_drift(grav, c->b, nuF, oi), gravity_drift(grav, c->b, nue, -oe)),
                 qe * 1e12);
    V3 jE = vmul(vsub(gravity_drift(grav, c->b, nuE, om), gravity_drift(grav, c->b, nue, -oe)),
                 qe * c->nE);
    c->cF[0] = vdot(jF, c->gx);
    c->cF[1] = vdot(jF, c->gy);
    c->cE[0] = vdot(jE, c->gx);
    c->cE[1] = vdot(jE, c->gy);
}
/**
 * Populate the explicitly idealized analytic height-dependent collision/conductivity background.
 */
void tube_conductivity(TubeGrid *g, double e_peak) {
    int n = g->c.nx * g->c.ny * g->c.ns;
    for (int i = 0; i < n; i++) {
        TubeCell *c = g->cell + i;
        double h = c->height;
        double nuF = 26.7 * exp(-h / 81400.) + 400 * exp(-(h - 110e3) / 10e3);
        double nue = 5e4 * exp(-(h - 110e3) / 20e3) + .01;
        conductivity_cell(c, e_peak, nuF, nue);
    }
}
/**
 * Read a strictly increasing height/collision table, require full geometry coverage, log-
 * interpolate collision rates and construct responses.
 */
int tube_conductivity_table(TubeGrid *g, double e_peak, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f)
        return -1;
    double height[4096], ion[4096], electron[4096];
    int n = 0;
    char line[512];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n' || !strncmp(line, "height_m,nu_O_s1,nu_e_s1", 22))
            continue;
        double h, a, b;
        char tail;
        if (n == 4096 || sscanf(line, "%lf,%lf,%lf %c", &h, &a, &b, &tail) != 3 || !isfinite(h) ||
            !isfinite(a) || !isfinite(b) || a <= 0 || b <= 0 || (n && h <= height[n - 1])) {
            fclose(f);
            return -2;
        }
        height[n] = h;
        ion[n] = log(a);
        electron[n] = log(b);
        n++;
    }
    fclose(f);
    if (n < 2)
        return -2;
    int cells = g->c.nx * g->c.ny * g->c.ns;
    for (int i = 0; i < cells; i++)
        if (g->cell[i].height < height[0] || g->cell[i].height > height[n - 1])
            return -3;
    for (int i = 0; i < cells; i++) {
        TubeCell *c = g->cell + i;
        int lo = 0, hi = n - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (height[mid] > c->height)
                hi = mid;
            else
                lo = mid;
        }
        double w = (c->height - height[lo]) / (height[hi] - height[lo]);
        conductivity_cell(c, e_peak, exp(ion[lo] + w * (ion[hi] - ion[lo])),
                          exp(electron[lo] + w * (electron[hi] - electron[lo])));
    }
    return 0;
}
/* Same response kernel as the legacy height table. The new input supplies
 * independently evaluated collision frequencies at every canonical cell. */
/**
 * Construct the same conductivity response from externally supplied per-canonical-cell collision
 * frequencies.
 */
void tube_conductivity_values(TubeGrid *g, double e_peak, const double *values) {
    size_t n = (size_t)g->c.nx * g->c.ny * g->c.ns;
    for (size_t i = 0; i < n; i++)
        conductivity_cell(g->cell + i, e_peak, values[6 * i + 4], values[6 * i + 5]);
}
