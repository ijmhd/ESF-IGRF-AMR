#ifndef ESF_SOLVER_H
#define ESF_SOLVER_H
#include <stddef.h>
#include <stdint.h>
#include "potential.h"
#include "electric_driver.h"
#include "integrator.h"
#ifdef ESF_MPI
#include <mpi.h>
#endif
typedef struct {
    GridConfig grid;
    int degree, parallel, hall, max_steps, regional_x;
    int petsc_potential, aligned_dipole, background_profile, periodic_x, open_apex, open_ends,
        parallel_order, progress_steps, parallel_relaxation;
    int amr_mode, amr_case, amr_coo, amr_levels, geometry_only;
    int amr_structure, amr_slab, amr_seed_guard, amr_fixture, amr_residual;
    double amr_curve_on, amr_curve_off, amr_velocity_on, amr_sigma_on, amr_quiet_s;
    double amr_seed_lambda, amr_branch_lambda, amr_seed_amplitude;
    double amr_work_weight;
    int amr_buffer, amr_hold;
    double amr_interval, amr_deplete_on, amr_deplete_off, amr_grad_on, amr_grad_off;
    int exb_pdm4;
    double pdm_beta;
    double dipole_equator_nT;
    double year, end, output, amplitude, seed_height, seed_width, e_peak, temperature, cfl, dtmax,
        rtol, background_east, seed_east;
    double seed_offset_deg; /* Apex-longitude offset from grid center; zero preserves prior
                               behavior. */
    int night_chemistry, night_reactions;
    int spatial_neutral;
    uint64_t neutral_grid_hash;
    double utc_start, chem_Ti, chem_Te;
    const ElectricDriver *driver;
    char driver_path[1024];
    char neutral_profile[1024], night_restart[1024];
    char coeff[1024], out[1024], collision_profile[1024], restart[1024], geometry_cache[1024];
} Config3;
enum {
    PROF_STATE,
    PROF_COEFF,
    PROF_COEFF_MPI,
    PROF_ASSEMBLY,
    PROF_SOLVE,
    PROF_BCAST,
    PROF_FACES,
    PROF_CELLS,
    PROF_REDUCE,
    PROF_VALID,
    PROF_COUNT
};
typedef struct {
    int peer, recv_first, recv_count, send_first, send_count;
    double *recv, *send;
} HaloPeer3;
typedef struct {
    Config3 c;
    TubeGrid g;
    Potential p;
    Igrf field;
    int rank, size, first, count, n, total, *counts, *offsets;
    double *rho, *momentum, *stage_n, *stage_m, *old_n, *old_m, *rhs_n, *rhs_m, *first_rhs_m;
    double *coeff, *local_coeff, *phi, *flux_n[3], *flux_m[3], *rate_face[3];
    double *cached_velocity, *cached_logrho, *pdm_weights[3];
    double cfl_rate, mass0, max_potential_residual, max_up, max_parallel, max_drift, wall0;
    double background_potential_gradient, boundary_rate, boundary_total, east_scale_m;
#ifdef ESF_MPI
    MPI_Request *halo_requests;
    MPI_Status *halo_statuses;
#endif
    void *petsc_context;
    HaloPeer3 *halo;
    int halo_peers, *coefficient_counts, *coefficient_offsets;
    double timing[PROF_COUNT];
    long potential_iterations;
    int solves, retries;
    V3 target_apex;
} Solver3;
_Noreturn void fail3(Solver3 *s, const char *message);
int solver3_load_field(const Config3 *c, Igrf *field);
double solver3_initial_density(const Config3 *c, double height);
void *alloc3(Solver3 *s, size_t count, size_t width);
void solver3_create(Solver3 *s, Config3 c, int rank, int size);
void solver3_free(Solver3 *s);
void solver3_gather(Solver3 *s, double *x);
int solver3_rhs(Solver3 *s, double *n, double *m, double *dn, double *dm);
int solver3_transport(Solver3 *s, const double *n, const double *m, double *dn, double *dm);
int solver3_valid(Solver3 *s, const double *n, const double *m);
double solver3_mass(Solver3 *s, const double *n);
void solver3_output(Solver3 *s, double time, int frame, int step);
double wall3(void);
void solver3_restart(Solver3 *s, const char *path, double *time, int *step);
uint64_t solver3_file_hash(Solver3 *s, const char *path);
void solver3_checkpoint_parameters(const Config3 *c, double expected[32]);
#ifdef ESF_PETSC
int amr3_run(Config3 c, int rank, int size);
int solver3_petsc_create(Solver3 *s);
int solver3_petsc_solve(Solver3 *s);
int solver3_petsc_seed(Solver3 *s);
void solver3_petsc_destroy(Solver3 *s);
#endif
#endif
