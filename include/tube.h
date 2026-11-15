#ifndef ESF_TUBE_H
#define ESF_TUBE_H
#include "igrf.h"
typedef struct {
    /* nx,ny: transverse logical cells; ns: whole conjugate tube segments.
     * lon/span are apex geographic longitude in degrees. All h and step
     * values are metres; h is spherical altitude above radius 6371200 m. */
    int nx, ny, ns;
    double lon, span, hmin, hmax, hbase, trace_step;
} GridConfig;
typedef struct {
    /* ECEF centre and a geometric field-aligned unit tangent. gx,gy are
     * contravariant gradients of normalized transverse labels, in 1/metre.
     * volume [m^3], B [T], height [m], gravity_parallel [m/s^2], nu [1/s]. */
    V3 position, b, gx, gy;
    double volume, B, height, gravity_parallel, nu;
    /* Conductivity and gravity-current responses per normalized F density,
     * plus a separately specified fixed molecular E-region population. */
    double pF, hF, pE, hE, cF[2], cE[2], nE;
} TubeCell;
typedef struct {
    /* e dot (d_phi/d_xi,d_phi/d_eta) is an oriented integrated E-cross-B
     * volume flux [m^3/s]. parallel_area [m^2] multiplies n*u_parallel;
     * magnetic_flux [T m^2] is a geometry-quality diagnostic. */
    double e[2], parallel_area, magnetic_flux;
} TubeFace;
typedef struct {
    GridConfig c;
    V3 *node;
    TubeCell *cell;
    TubeFace *face[3];
    double max_alignment_error, max_parallel_electric_metric_error, min_volume,
        max_magnetic_flux_jump;
} TubeGrid;
int tube_geometry(TubeGrid *g, const Igrf *m, GridConfig c);
int tube_geometry_allocate(TubeGrid *g, GridConfig c);
int tube_geometry_rows(TubeGrid *g, const Igrf *m, int first, int last);
int tube_cell_moments(const TubeGrid *g, int i, int j, int k, double out[12]);
int tube_moment_weights(const double moments[4][4], double weights[4]);
int tube_pdm_weights(const TubeGrid *g, int periodic_x, double *weights[3]);
double tube_pdm_limit(double high, double upstream, double donor, double downstream, double beta);
void tube_free(TubeGrid *g);
int tube_apex(const Igrf *m, V3 start, V3 *apex);
int tube_apex_label(const Igrf *m, double lon, double height, V3 *apex);
int tube_trace(const Igrf *m, double lon, double height, double base, double step, int ns, V3 *out);
int tube_cell_id(const TubeGrid *g, int i, int j, int k);
int tube_face_id(const TubeGrid *g, int axis, int i, int j, int k);
int tube_face_count(const TubeGrid *g, int axis);
double chapman(double height, double peak, double scale);
void tube_conductivity(TubeGrid *g, double e_peak);
int tube_conductivity_table(TubeGrid *g, double e_peak, const char *path);
#endif
