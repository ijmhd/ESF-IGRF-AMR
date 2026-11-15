/*
 * Small IGRF command-line diagnostic; angles are degrees, positions metres and magnetic output
 * nT.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "igrf.h"
#include <stdio.h>
#include <stdlib.h>
/**
 * Parse this executable's interface, validate inputs, run its diagnostic/model entry point, and
 * return a nonzero status on failure.
 */
int main(int argc, char **argv) {
    if (argc != 7) {
        fprintf(stderr, "field_eval COEFF YEAR LAT LON ALT_KM DEGREE\n");
        return 2;
    }
    Igrf m;
    if (igrf_load(&m, argv[1], atof(argv[2]), atoi(argv[6])))
        return 1;
    double lat = atof(argv[3]), lon = atof(argv[4]), h = atof(argv[5]) * 1000;
    V3 p = geodetic_position(lat, lon, h), b = igrf_field(&m, p),
       ned = vmul(field_ned(b, lat, lon), 1e9);
    printf("%.15g %.15g %.15g %.15g\n", ned.x, ned.y, ned.z, vnorm(ned));
    return 0;
}
