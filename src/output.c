/*
 * Uniform-grid metadata, field outputs and checkpoints. Always use saved metadata to decode
 * shapes and units.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/**
 * Open a named uniform-grid output file in the run directory; abort on an I/O failure.
 */
static FILE *open3(Solver3 *s, const char *name, const char *mode) {
    char path[1200];
    snprintf(path, sizeof(path), "%s/%s", s->c.out, name);
    FILE *f = fopen(path, mode);
    if (!f)
        fail3(s, "cannot write output");
    return f;
}
/**
 * Write exactly n binary double values and reject short writes.
 */
static void doubles(FILE *f, const double *x, size_t n) {
    unsigned short one = 1;
    if (*(unsigned char *)&one) {
        if (fwrite(x, sizeof(double), n, f) != n)
            abort();
    } else
        for (size_t i = 0; i < n; i++) {
            unsigned char b[8], r[8];
            memcpy(b, x + i, 8);
            for (int j = 0; j < 8; j++)
                r[j] = b[7 - j];
            if (fwrite(r, 8, 1, f) != 1)
                abort();
        }
}
/**
 * Copy/hash the effective collision input so the saved run records its actual background policy.
 */
static void collision_input(Solver3 *s) {
    if (!s->c.collision_profile[0])
        return;
    FILE *in = fopen(s->c.collision_profile, "rb"), *out = open3(s, "collision_profile.csv", "wb");
    if (!in)
        fail3(s, "cannot archive collision input");
    unsigned char block[4096];
    size_t n;
    while ((n = fread(block, 1, sizeof(block), in)))
        if (fwrite(block, 1, n, out) != n)
            fail3(s, "cannot archive collision input");
    if (ferror(in))
        fail3(s, "cannot read collision input");
    fclose(in);
    fclose(out);
    out = open3(s, "collision_input.json", "w");
    fputs(
        "{\"file\":\"collision_profile.csv\",\"format\":\"height_m,nu_O_s1,nu_e_s1\",\"interpolation\":\"linear in log frequency\",\"scope\":\"same O+ collision frequency used in parallel drag, conductivity and gravity current; constant fluid temperature remains specified separately\"}\n",
        out);
    fclose(out);
}
/**
 * Gather uniform conserved states, write fields/metadata/diagnostics and an atomically replaced
 * portable checkpoint.
 */
void solver3_output(Solver3 *s, double time, int frame, int step) {
#ifdef ESF_MPI
    MPI_Gatherv(s->rank ? s->rho + s->first * s->c.grid.ns : MPI_IN_PLACE, s->n, MPI_DOUBLE, s->rho,
                s->counts, s->offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(s->rank ? s->momentum + s->first * s->c.grid.ns : MPI_IN_PLACE, s->n, MPI_DOUBLE,
                s->momentum, s->counts, s->offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    double mass = solver3_mass(s, s->rho);
    if (s->rank)
        return;
    int nx = s->c.grid.nx, ny = s->c.grid.ny, ns = s->c.grid.ns;
    char name[128];
    if (!frame) {
        FILE *f = open3(s, "metadata.json", "w");
        fprintf(
            f,
            "{\n\"model\":\"independent field-aligned 3D electrostatic ESF research prototype\",\n"
            "\"field\":\"%s\",\"igrf_generation\":%d,\"degree\":%d,\"year\":%.17g,\"ranks\":%d,\n"
            "\"nx\":%d,\"ny\":%d,\"ns\":%d,\"longitude_center_deg\":%.17g,\"longitude_span_deg\":%.17g,\n"
            "\"apex_min_m\":%.17g,\"apex_max_m\":%.17g,\"foot_spherical_altitude_m\":%.17g,\"reference_radius_m\":6371200,\n"
            "\"trace_step_max_m\":%.17g,\"parallel_flow\":%d,\"hall\":%d,\"Ti_plus_Te_K\":%.17g,\n"
            "\"initial_seed_amplitude\":%.17g,\"seed_apex_height_m\":%.17g,\"seed_apex_width_m\":%.17g,\"fixed_E_peak_density_m3\":%.17g,\n"
            "\"seed_east_width_m\":%.17g,\"seed_longitudinal_profile\":\"cosine when width=0; (2u^2-1)exp(-u^2) localized packet otherwise\",\n"
            "\"cfl\":%.17g,\"dt_max_s\":%.17g,\"rtol\":%.17g,\"petsc_potential\":%d,\n"
            "\"background_east_apex_V_m\":%.17g,\"background_potential_gradient_V\":%.17g,\n"
            "\"regional_x\":%d,\"regional_boundary_model_version\":2,\n"
            "\"periodic_x\":%d,\"open_apex\":%d,\"open_field_line_ends\":%d,\"parallel_reconstruction_order\":%d,\n"
            "\"initial_profile\":\"%s\",\"ideal_dipole_equator_nT\":%.17g,\n"
            "\"parallel_mass_flux\":\"%s\",\n"
            "\"parallel_relaxation_momentum\":\"mass-consistent velocity advection and logarithmic isothermal pressure source; not a shock-conservative Euler solver\",\n"
            "\"background_E_calibration\":\"exact seed-height and central-label apex dual basis, independent of mesh centers\",\n"
            "\"ExB_reconstruction\":\"%s\",\"PDM_beta\":%.17g,\n"
            "\"PDM_boundary_treatment\":\"MUSCL in two physical boundary layers; periodic longitude uses wrapped cubic moments\",\n"
            "\"transport\":\"conservative 3D FV; configurable volume-moment PDM4 or MUSCL ExB transport, parallel isothermal Rusanov with configurable reconstruction order; ETD2 (SSPRK2 density) with stiff ion-neutral drag\",\n"
            "\"potential\":\"Q1 field-line-integrated variable Pedersen/Hall tensor; native ILU0/BiCGStab or distributed PETSc/FGMRES selected by petsc_potential; true residual checked\",\n"
            "\"boundaries\":\"zero polarization potential at apex-height edges; longitude periodic if periodic_x, otherwise regional-x uses d_xi(phi_pol)=0 including zero imposed E, otherwise zero polarization; regional longitude particles use zero-order open extrapolation; remaining particles follow explicit open flags or imposed-E apex settings\",\n"
            "\"parallel_electric_field\":\"electrostatic perturbation equipotential along traced lines; ambipolar pressure represented in parallel fluid equation\",\n"
            "\"parallel_inertia_model\":\"advected parallel scalar velocity; perpendicular curvature inertia neglected\",\n"
            "\"density_floor\":false,\"conductivity_floor\":false,\"velocity_cap\":false,\"amr\":false,\n"
            "\"background_status\":\"Chapman F fluid; constant fluid temperature; analytical collision background unless collision_input.json supplies an external profile; no event-specific calibration\",\n"
            "\"mpi_strategy\":\"owned full flux-tube columns; neighbor state halos; owned coefficients gathered (allgathered for PETSc); native root solve or PETSc distributed rows; full state gathered only for output; geometry replicated\",\n"
            "\"density_unit_m3\":1e12,\"geometry_cell_columns\":[\"ecef_x_m\",\"ecef_y_m\",\"ecef_z_m\",\"spherical_altitude_m\",\"B_T\",\"volume_m3\"],\n"
            "\"binary_format\":\"little-endian float64; geometry_cells [ny,nx,ns,6]; nodes [ny+1,nx+1,ns+1,3]; states [2,ny,nx,ns] = rho and rho*u_parallel; phi [ny+1,nx+1] volts\",\n"
            "\"geometry_quality\":{\"field_direction_sine_error_max\":%.17g,\"E_parallel_metric_relative_max\":%.17g,\"adjacent_parallel_face_magnetic_flux_relative_jump_max\":%.17g,\"minimum_volume_m3\":%.17g},\n"
            "\"hong_kong_target\":{\"geodetic_lat_deg\":22.3,\"lon_deg\":114.2,\"altitude_m\":350000,\"field_line_apex_spherical_altitude_m\":%.17g,\"field_line_apex_lon_deg\":%.17g}\n}\n",
            s->c.aligned_dipole ? "axis-aligned ideal dipole"
                                : (s->c.degree == 1 ? "tilted degree-1 dipole" : "IGRF-14"),
            s->c.aligned_dipole ? 0 : 14, s->c.degree, s->c.year, s->size, nx, ny, ns,
            s->c.grid.lon, s->c.grid.span, s->c.grid.hmin, s->c.grid.hmax, s->c.grid.hbase,
            s->c.grid.trace_step, s->c.parallel, s->c.hall, s->c.temperature, s->c.amplitude,
            s->c.seed_height, s->c.seed_width, s->c.e_peak, s->c.seed_east, s->c.cfl, s->c.dtmax,
            s->c.rtol, s->c.petsc_potential, s->c.background_east, s->background_potential_gradient,
            s->c.regional_x, s->c.periodic_x, s->c.open_apex, s->c.open_ends, s->c.parallel_order,
            s->c.background_profile ? "TELF double Chapman O+" : "single Chapman O+",
            s->c.dipole_equator_nT,
            s->c.parallel_relaxation
                ? "quadratic relaxation blend of Rusanov and Scharfetter-Gummel"
                : "Rusanov",
            s->c.exb_pdm4 ? "PDM4 with trilinear physical-volume moments" : "MUSCL", s->c.pdm_beta,
            s->g.max_alignment_error, s->g.max_parallel_electric_metric_error,
            s->g.max_magnetic_flux_jump, s->g.min_volume, vnorm(s->target_apex) - 6371200.,
            atan2(s->target_apex.y, s->target_apex.x) * 180 / 3.14159265358979323846);
        fclose(f);
        collision_input(s);
        f = open3(s, "geometry_cells.bin", "wb");
        for (int i = 0; i < s->total; i++) {
            TubeCell *c = s->g.cell + i;
            double v[6] = {c->position.x, c->position.y, c->position.z, c->height, c->B, c->volume};
            doubles(f, v, 6);
        }
        fclose(f);
        f = open3(s, "nodes.bin", "wb");
        for (int i = 0; i < (nx + 1) * (ny + 1) * (ns + 1); i++) {
            V3 x = s->g.node[i];
            double v[3] = {x.x, x.y, x.z};
            doubles(f, v, 3);
        }
        fclose(f);
        f = open3(s, "diagnostics.csv", "w");
        fputs(
            "step,time_s,min_rho,max_rho,mass,relative_mass_error,max_up_m_s,max_ExB_m_s,max_parallel_m_s,potential_residual,potential_iterations,solves,retries,wall_s,integrated_boundary_particle_flux\n",
            f);
        fclose(f);
        f = open3(s, "frames.csv", "w");
        fputs("frame,time_s,state,phi,columns\n", f);
        fclose(f);
    }
    snprintf(name, sizeof(name), "state_%04d.bin", frame);
    FILE *f = open3(s, name, "wb");
    doubles(f, s->rho, s->total);
    doubles(f, s->momentum, s->total);
    fclose(f);
    snprintf(name, sizeof(name), "phi_%04d.bin", frame);
    f = open3(s, name, "wb");
    doubles(f, s->phi, (nx + 1) * (ny + 1));
    fclose(f);
    snprintf(name, sizeof(name), "columns_%04d.csv", frame);
    f = open3(s, name, "w");
    fputs(
        "i,j,apex_label_m,longitude_label_deg,apex_rho,pedersen_F_S,pedersen_E_S,hall_F_S,hall_E_S,HK_latitude_cut_lon_deg,HK_latitude_cut_altitude_m,HK_latitude_cut_rho\n",
        f);
    double low = HUGE_VAL, high = 0, up = -HUGE_VAL, drift = 0, parallel = 0;
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            double d0 = .5 * nx *
                        (s->phi[j * (nx + 1) + i + 1] + s->phi[(j + 1) * (nx + 1) + i + 1] -
                         s->phi[j * (nx + 1) + i] - s->phi[(j + 1) * (nx + 1) + i]);
            double d1 = .5 * ny *
                        (s->phi[(j + 1) * (nx + 1) + i] + s->phi[(j + 1) * (nx + 1) + i + 1] -
                         s->phi[j * (nx + 1) + i] - s->phi[j * (nx + 1) + i + 1]);
            double sigma[4] = {0, 0, 0, 0}, hmax = -HUGE_VAL, apex_n = 0, cut_h = NAN,
                   cut_lon = NAN, cut_n = NAN;
            for (int k = 0; k < ns; k++) {
                int id = tube_cell_id(&s->g, i, j, k);
                TubeCell *c = s->g.cell + id;
                double n = s->rho[id];
                V3 E = vmul(vadd(vmul(c->gx, d0), vmul(c->gy, d1)), -1),
                   v = vmul(vcross(E, c->b), 1 / c->B);
                up = fmax(up, vdot(v, vunit(c->position)));
                drift = fmax(drift, vnorm(v));
                parallel = fmax(parallel, n > 0 ? fabs(s->momentum[id] / n) : 0);
                low = fmin(low, n);
                high = fmax(high, n);
                double area =
                           .5 * (s->g.face[2][tube_face_id(&s->g, 2, i, j, k)].parallel_area +
                                 s->g.face[2][tube_face_id(&s->g, 2, i, j, k + 1)].parallel_area),
                       ds = c->volume / area;
                sigma[0] += c->pF * n * ds;
                sigma[1] += c->pE * ds;
                sigma[2] += c->hF * n * ds;
                sigma[3] += c->hE * ds;
                if (c->height > hmax) {
                    hmax = c->height;
                    apex_n = n;
                }
                if (k) {
                    TubeCell *prev = c - 1;
                    double la0, lo0, ha0, la1, lo1, ha1;
                    geodetic_coordinates(prev->position, &la0, &lo0, &ha0);
                    geodetic_coordinates(c->position, &la1, &lo1, &ha1);
                    if (la0 <= 22.3 && la1 >= 22.3) {
                        double w = (22.3 - la0) / (la1 - la0);
                        cut_h = ha0 + w * (ha1 - ha0);
                        cut_lon = lo0 + w * (lo1 - lo0);
                        cut_n = s->rho[id - 1] + w * (n - s->rho[id - 1]);
                    }
                }
            }
            fprintf(f, "%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n", i, j,
                    s->c.grid.hmin + (s->c.grid.hmax - s->c.grid.hmin) * (j + .5) / ny,
                    s->c.grid.lon + s->c.grid.span * ((i + .5) / nx - .5), apex_n, sigma[0],
                    sigma[1], sigma[2], sigma[3], cut_lon, cut_h, cut_n);
        }
    fclose(f);
    s->max_up = up;
    s->max_parallel = parallel;
    s->max_drift = drift;
    f = open3(s, "diagnostics.csv", "a");
    fprintf(f, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%.6f,%.17g\n",
            step, time, low, high, mass, (mass - s->mass0 - s->boundary_total) / s->mass0, up,
            drift, parallel, s->p.residual, s->p.iterations, s->solves, s->retries,
            wall3() - s->wall0, s->boundary_total);
    fclose(f);
    f = open3(s, "frames.csv", "a");
    fprintf(f, "%d,%.17g,state_%04d.bin,phi_%04d.bin,columns_%04d.csv\n", frame, time, frame, frame,
            frame);
    fclose(f);
    printf("t=%7.1f cells=%d min=%.3g up=%.3g parallel=%.3g residual=%.2g mass=%.2g wall=%.1fs\n",
           time, s->total, low, up, parallel, s->p.residual,
           (mass - s->mass0 - s->boundary_total) / s->mass0, wall3() - s->wall0);
    fflush(stdout);
}
