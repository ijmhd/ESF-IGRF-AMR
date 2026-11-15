/*
 * Geometry inspection program; it traces tubes and reports the Hong Kong target mapping,
 * independently of plasma evolution.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* Geodetic target -> apex-labelled complete conjugate flux tube.
 * Geometry diagnostic only; no density evolution or magnetic coordinates are
 * inferred from a dipole approximation. */
#include "tube.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Parse this executable's interface, validate inputs, run its diagnostic/model entry point, and
 * return a nonzero status on failure.
 */
int main(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr, "tube_map COEFF YEAR LAT_DEG LON_DEG ALT_KM DEGREE LINE_CSV\n");
        return 2;
    }
    char *end;
    double value[5];
    for (int i = 0; i < 5; i++) {
        value[i] = strtod(argv[i + 2], &end);
        if (*end || end == argv[i + 2] || !isfinite(value[i]))
            return 2;
    }
    if (value[1] < -89 || value[1] > 89 || value[3] < 80 || value[3] > 3000 ||
        value[4] != floor(value[4]))
        return 2;
    Igrf m;
    if (igrf_load(&m, argv[1], value[0], (int)value[4]))
        return 1;
    V3 target = geodetic_position(value[1], value[2], value[3] * 1000), apex;
    if (tube_apex(&m, target, &apex))
        return 1;
    double h = vnorm(apex) - 6371200., lon = atan2(apex.y, apex.x) * 180 / 3.14159265358979323846;
    enum { N = 1024 };
    V3 line[N + 1];
    if (tube_trace(&m, lon, h, 90000, 1000, N, line))
        return 1;
    FILE *f = fopen(argv[7], "w");
    if (!f)
        return 1;
    fputs(
        "k,ecef_x_m,ecef_y_m,ecef_z_m,geodetic_lat_deg,lon_deg,geodetic_alt_m,spherical_alt_m,B_T\n",
        f);
    double distance = HUGE_VAL;
    for (int i = 0; i <= N; i++) {
        double la, lo, alt;
        geodetic_coordinates(line[i], &la, &lo, &alt);
        fprintf(f, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n", i, line[i].x, line[i].y,
                line[i].z, la, lo, alt, vnorm(line[i]) - 6371200., vnorm(igrf_field(&m, line[i])));
        if (i) {
            V3 d = vsub(line[i], line[i - 1]);
            double t = vdot(vsub(target, line[i - 1]), d) / vdot(d, d);
            t = fmax(0, fmin(1, t));
            distance = fmin(distance, vnorm(vsub(target, vadd(line[i - 1], vmul(d, t)))));
        }
    }
    fclose(f);
    V3 b = vmul(field_ned(igrf_field(&m, target), value[1], value[2]), 1e9);
    printf(
        "{\"year\":%.17g,\"degree\":%d,\"target_geodetic_lat_deg\":%.17g,\"target_lon_deg\":%.17g,\"target_altitude_m\":%.17g,\"apex_spherical_altitude_m\":%.17g,\"apex_lon_deg\":%.17g,\"apex_geocentric_lat_deg\":%.17g,\"target_NED_nT\":[%.17g,%.17g,%.17g],\"target_distance_to_traced_polyline_m\":%.17g,\"base_spherical_altitude_m\":90000}\n",
        value[0], m.degree, value[1], value[2], value[3] * 1000, h, lon,
        atan2(apex.z, hypot(apex.x, apex.y)) * 180 / 3.14159265358979323846, b.x, b.y, b.z,
        distance);
    return 0;
}
