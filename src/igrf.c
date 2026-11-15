/*
 * Independent Schmidt semi-normalized IGRF evaluator. Position is ECEF metres, returned field is
 * ECEF tesla. Geodetic WGS84 height differs from spherical radius minus 6371.2 km.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
/* Independently derived Schmidt semi-normalized spherical harmonics.
 * Coefficients: official IAGA/NOAA IGRF-14; field output is tesla, ECEF.
 * Reference radius 6371.2 km differs from the WGS84 geodetic ellipsoid. */
#include "igrf.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.14159265358979323846
/**
 * Return component-wise ECEF vector addition.
 */
V3 vadd(V3 a, V3 b) {
    return (V3){a.x + b.x, a.y + b.y, a.z + b.z};
}
/**
 * Return component-wise ECEF vector subtraction.
 */
V3 vsub(V3 a, V3 b) {
    return (V3){a.x - b.x, a.y - b.y, a.z - b.z};
}
/**
 * Multiply an ECEF vector by a scalar.
 */
V3 vmul(V3 a, double s) {
    return (V3){a.x * s, a.y * s, a.z * s};
}
/**
 * Return the Euclidean dot product.
 */
double vdot(V3 a, V3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
/**
 * Return the right-handed Euclidean cross product.
 */
V3 vcross(V3 a, V3 b) {
    return (V3){a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
/**
 * Return the Euclidean vector magnitude.
 */
double vnorm(V3 a) {
    return sqrt(vdot(a, a));
}
/**
 * Normalize a nonzero vector; return the zero vector for zero magnitude.
 */
V3 vunit(V3 a) {
    double n = vnorm(a);
    return n ? vmul(a, 1 / n) : (V3){0, 0, 0};
}
/**
 * Read all 195 official IGRF-14 Gauss coefficients, interpolate five-year epochs or apply
 * 2025-2030 secular variation, and retain the requested harmonic degree.
 */
int igrf_load(Igrf *m, const char *file, double year, int degree) {
    if (year < 1900 || year > 2030 || degree < 1 || degree > 13)
        return -1;
    FILE *fp = fopen(file, "r");
    if (!fp)
        return -2;
    memset(m, 0, sizeof(*m));
    m->year = year;
    m->degree = degree;
    char line[4096];
    int count = 0;
    while (fgets(line, sizeof(line), fp)) {
        if ((line[0] != 'g' && line[0] != 'h') || line[1] != ' ')
            continue;
        char *t = strtok(line, " \t\r\n");
        char kind = *t;
        t = strtok(NULL, " \t\r\n");
        if (!t)
            continue;
        int n = atoi(t);
        t = strtok(NULL, " \t\r\n");
        if (!t)
            continue;
        int k = atoi(t);
        if (n < 1 || n > 13 || k < 0 || k > n)
            continue;
        double values[27];
        int nv = 0;
        while ((t = strtok(NULL, " \t\r\n")) && nv < 27) {
            char *end;
            double v = strtod(t, &end);
            if (*end)
                break;
            values[nv++] = v;
        }
        if (nv != 27) {
            fclose(fp);
            return -3;
        }
        int i = (int)floor((fmin(year, 2025.) - 1900.) / 5.);
        double value = i < 25
                           ? values[i] + (year - 1900. - 5 * i) / 5. * (values[i + 1] - values[i])
                           : values[25] + (year - 2025.) * values[26];
        if (kind == 'g')
            m->g[n][k] = value;
        else
            m->h[n][k] = value;
        count++;
    }
    fclose(fp);
    return count == 195 ? 0 : -4;
}
/**
 * Evaluate Schmidt semi-normalized harmonics and their derivatives, convert spherical field to
 * ECEF tesla, and use a continuous near-pole evaluation.
 */
V3 igrf_field(const Igrf *m, V3 position) {
    double r = vnorm(position), theta = atan2(hypot(position.x, position.y), position.z);
    double lon = atan2(position.y, position.x), s = sin(theta), c = cos(theta);
    /* Stable evaluation at a pole using its continuous Cartesian limit. */
    double evaltheta = theta;
    if (fabs(s) < 1e-10) {
        evaltheta = theta < PI / 2 ? 1e-10 : PI - 1e-10;
        s = sin(evaltheta);
        c = cos(evaltheta);
    }
    double p[14][14] = {{0}}, d[14][14] = {{0}};
    p[0][0] = 1.;
    for (int k = 0; k <= m->degree; k++) {
        if (k) {
            p[k][k] = (2 * k - 1) * s * p[k - 1][k - 1];
            d[k][k] = (2 * k - 1) * (c * p[k - 1][k - 1] + s * d[k - 1][k - 1]);
        }
        if (k < m->degree) {
            p[k + 1][k] = (2 * k + 1) * c * p[k][k];
            d[k + 1][k] = (2 * k + 1) * (-s * p[k][k] + c * d[k][k]);
        }
        for (int n = k + 2; n <= m->degree; n++) {
            p[n][k] = ((2 * n - 1) * c * p[n - 1][k] - (n + k - 1) * p[n - 2][k]) / (n - k);
            d[n][k] =
                ((2 * n - 1) * (-s * p[n - 1][k] + c * d[n - 1][k]) - (n + k - 1) * d[n - 2][k]) /
                (n - k);
        }
    }
    double br = 0, bt = 0, bp = 0;
    for (int n = 1; n <= m->degree; n++) {
        double radial = pow(6371200. / r, n + 2);
        for (int k = 0; k <= n; k++) {
            double norm = 1.;
            if (k) {
                norm = 2.;
                for (int j = n - k + 1; j <= n + k; j++)
                    norm /= j;
                norm = sqrt(norm);
            }
            double cl = cos(k * lon), sl = sin(k * lon), g = m->g[n][k], h = m->h[n][k];
            br += (n + 1) * radial * (g * cl + h * sl) * p[n][k] * norm;
            bt -= radial * (g * cl + h * sl) * d[n][k] * norm;
            bp += radial * k * (g * sl - h * cl) * p[n][k] * norm / s;
        }
    }
    V3 er = {sin(theta) * cos(lon), sin(theta) * sin(lon), cos(theta)};
    V3 et = {cos(theta) * cos(lon), cos(theta) * sin(lon), -sin(theta)},
       ep = {-sin(lon), cos(lon), 0};
    return vmul(vadd(vadd(vmul(er, br), vmul(et, bt)), vmul(ep, bp)), 1e-9);
}
/**
 * Convert geocentric latitude/longitude in degrees and radius in metres to ECEF metres.
 */
V3 spherical_position(double lat, double lon, double r) {
    lat *= PI / 180.;
    lon *= PI / 180.;
    return (V3){r * cos(lat) * cos(lon), r * cos(lat) * sin(lon), r * sin(lat)};
}
/**
 * Convert WGS84 geodetic latitude/longitude/altitude to ECEF metres.
 */
V3 geodetic_position(double lat, double lon, double h) {
    lat *= PI / 180.;
    lon *= PI / 180.;
    double e2 = 6.6943799901413165e-3;
    double n = 6378137. / sqrt(1 - e2 * sin(lat) * sin(lat));
    return (V3){(n + h) * cos(lat) * cos(lon), (n + h) * cos(lat) * sin(lon),
                (n * (1 - e2) + h) * sin(lat)};
}
/**
 * Recover WGS84 geodetic latitude/longitude/altitude from ECEF coordinates by iterative latitude
 * refinement.
 */
void geodetic_coordinates(V3 p, double *lat, double *lon, double *h) {
    double e2 = 6.6943799901413165e-3, r = hypot(p.x, p.y), a = atan2(p.z, r * (1 - e2));
    for (int i = 0; i < 12; i++) {
        double n = 6378137. / sqrt(1 - e2 * sin(a) * sin(a));
        a = atan2(p.z + e2 * n * sin(a), r);
    }
    double n = 6378137. / sqrt(1 - e2 * sin(a) * sin(a));
    *h = r * cos(a) + p.z * sin(a) - n * (1 - e2 * sin(a) * sin(a));
    *lat = a * 180 / PI;
    *lon = atan2(p.y, p.x) * 180 / PI;
}
/**
 * Project an ECEF field into local north/east/down components at specified geographic angles.
 */
V3 field_ned(V3 b, double lat, double lon) {
    lat *= PI / 180;
    lon *= PI / 180;
    return (V3){vdot(b, (V3){-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat)}),
                vdot(b, (V3){-sin(lon), cos(lon), 0}),
                -vdot(b, (V3){cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat)})};
}
