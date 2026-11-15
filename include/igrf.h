#ifndef ESF_IGRF_H
#define ESF_IGRF_H
#define IGRF_N 13
typedef struct {
    double x, y, z;
} V3;
typedef struct {
    double g[14][14], h[14][14], year;
    int degree;
} Igrf;
V3 vadd(V3 a, V3 b);
V3 vsub(V3 a, V3 b);
V3 vmul(V3 a, double s);
double vdot(V3 a, V3 b);
V3 vcross(V3 a, V3 b);
double vnorm(V3 a);
V3 vunit(V3 a);
int igrf_load(Igrf *m, const char *file, double year, int degree);
V3 igrf_field(const Igrf *m, V3 ecef_m);
V3 spherical_position(double latitude_deg, double longitude_deg, double radius_m);
V3 geodetic_position(double latitude_deg, double longitude_deg, double altitude_m);
void geodetic_coordinates(V3 p, double *latitude_deg, double *longitude_deg, double *altitude_m);
V3 field_ned(V3 field, double latitude_deg, double longitude_deg);
#endif
