/* SPDX-License-Identifier: GPL-3.0-only */
#ifndef ESF_AMR3_H
#define ESF_AMR3_H
#include "solver.h"
#include <petscksp.h>
#include <stdio.h>

typedef struct {
    /* level is PRIVATE: -1/0/1 -> physical L0/L1/L2 for three levels.
     * i,j label this backing column; root is its containing internal root.
     * Actual density ownership also depends on k and the segmented masks. */
    int level, i, j, root;
} A3Leaf;
typedef struct {
    int axis, i, j, first, last, left, right, level, ll, rl;
} A3FaceBlock;
typedef struct {
    double volume, nu, gravity, bg, area[2];
} A3Metric;
typedef struct {
    /* Coefficients of a canonical potential value in independent DOFs:
     * phi_node = sum(w[q]*phi[id[q]]) + fixed + drive*external_E.
     * Hanging nodes therefore do not create independent unknowns. */
    int n, id[4];
    double w[4], fixed, drive;
} A3Map;
typedef struct {
    int peer, recv_first, recv_count, send_first, send_count;
    double *recv, *send;
} A3Halo;
typedef struct {
    /* The public grid dimensions are doubled internally for three levels.
     * Canonical geometry is then twice the INTERNAL nx,ny, giving four
     * times the user L0 resolution; ns is never refined. */
    Config3 c;
    TubeGrid *g, layout_storage[3], *layout;
    struct A3Night *night; /* Transported molecular species on canonical geometry. */
    int rank, size, nx, ny, fx, fy, ns, nb, nleaf, split;
    /* Owned rows are [row_first,row_last). Owned backing slots come first;
     * slots includes ghosts. global[local_slot] and slot[global_column]
     * translate the two index spaces; slot is -1 outside this rank's halo. */
    int row_first, row_last, owned_first, owned_count, slots, peers;
    /* refined controls the column union/field backing; density_refined and
     * macro_mask determine actual density DOFs at each along-field index.
     * ready contains lazy reconstruction validity bits, cleared each stage. */
    unsigned char *refined, *ready;
    int *cache_refs;
    int *first, *cover, *slot, *global, *row, *virtual_offset, *partition;
    A3Leaf *leaf;
    A3Metric *metric;
    A3Halo *halo;
    MPI_Request *requests;
    MPI_Status *statuses; /* Real storage avoids vendor-header STATUS_IGNORE diagnostics. */
    double *base_n, *base_m, *sx_n, *sy_n, *amp_n;
    double *root_volume, *root_bg, *fine_bg, *offset_x, *offset_y, *response;
    double *base_velocity, *velocity, *log_n, *child_n, *child_m, *child_velocity;
    double *weight_storage[3][3];
    double *(*weights)[3];
    /* Borrowed current-stage views, never independently freed. n is in
     * 1e12 m^-3 units; m=n*u_parallel [normalized density m/s]. */
    double *current_n, *current_m, field_min, field_max, material_min, material_max;
    A3Map *map;
    int ndof;
    double *phi, *rates, *transverse_rates, *parallel_n, *parallel_m, *parallel_rate;
    Mat matrix;
    Vec rhs, solution, residual, all;
    VecScatter scatter;
    KSP ksp;
    void *coo_context;
    FILE *audit_file;
    int solves, last_iterations;
    double potential_residual, max_residual;
    /* rate bounds explicit transport [1/s]; transverse_rate integrates the
     * grid-travel trigger. They are recomputed from actual current fluxes.
     * Neither is a physical plasma growth rate. */
    double rate, transverse_rate, remap_error, wall0, east_scale_m;
    int mesh_epoch, output_mesh_epoch, regrids, refinements, coarsenings, regrid_log_started;
    int mask_count, active_cells;
    unsigned char *density_refined, *seed_mask, *field_requested;
    /* Detector arrays describe mesh decisions only. In particular their
     * normalization scale is not a physical density/conductivity floor. */
    unsigned short *tag_reason;
    double *quiet_since, *born, *tag_values, *root_metric_store;
    double *marker, *base_marker, *marker_rhs;
    double tag_time, tag_start, tag_elapsed, seed_lambda, branch_lambda, density_scale;
    unsigned long unresolved_seed;
    V3 *tag_drift_fine, *tag_drift_root;
    double *tag_tensor;
    int tag_row_first, tag_row_last;
    double *parallel_length_storage[3], *parallel_response_storage[3];
    double **parallel_length, **parallel_response;
    unsigned char *macro_mask, *macro_next, *macro_field;
    double *macro_quiet, *macro_born, *macro_volume, *macro_bg, *macro_state[4];
    double *potential_gradient[3];
    double *etd_weight[3]; /* Owned-cell decay, phi1, phi2; valid for one attempted dt. */
    A3FaceBlock *face_blocks;
    int face_block_count, face_block_capacity;

    double time; /* actual SSPRK stage time, seconds from model UTC start */
    double t_regrid, max_potential_remap_change;
    double boundary_rate, boundary_total, import_mass_error, geometry_density_adjustment;
    /* Rank-local cumulative timers; regrid includes nested phases. They are
     * not mutually exclusive and must not be summed into a percent chart. */
    double t_exchange, t_prepare, t_coeff, t_potential, t_transport;
    unsigned long long rhs_calls, tag_checks, chemistry_calls, timestep_attempts, cfl_adjustments;
} A3;

_Noreturn void a3_fail(A3 *a, const char *message);
void *a3_alloc(A3 *a, size_t n, size_t width);
double a3_proposed_child(A3 *a, int b, int child, int k, int field);
int a3_owner(const A3 *a, int p, int k);
int a3_active(const A3 *a, int p, int k);
A3Leaf a3_density_leaf(const A3 *a, int p, int k);
const unsigned char *a3_density_mask(const A3 *a);
void a3_structure_tags(A3 *a, double *n, double *m, unsigned char *mask, int *quiet);
void a3_structure_output(A3 *a, const double *n, const double *marker, int frame, double time);
void a3_mesh(A3 *a, TubeGrid *g, Config3 c, int rank, int size, int mode);
void a3_mesh_mask(A3 *a, TubeGrid *g, Config3 c, int rank, int size, int mode,
                  const unsigned char *mask, A3 *reuse, const int *cuts);
void a3_mesh_mask_work(A3 *a, TubeGrid *g, Config3 c, int rank, int size, int mode,
                       const unsigned char *mask, A3 *reuse, const int *cuts,
                       const unsigned char *work_mask);
void a3_tags(A3 *a, double *n, double *m, unsigned char *mask, int *quiet);
void a3_free(A3 *a);
void a3_prepare(A3 *a, double *n, double *m);
double a3_sample(A3 *a, int field, int level, int i, int j, int k);
void a3_face_pair(A3 *a, int level, int axis, int i, int j, int k, int sign, double *n,
                  double *second);
void a3_initialize_marker(A3 *a, const double *n);
void a3_initialize(A3 *a, double *n, double *m, int test);
void a3_remap(A3 *old, A3 *next, double *n, double *m, double *new_n, double *new_m);
double a3_mass(A3 *a, const double *field);
int a3_valid(A3 *a, const double *n, const double *m);
int a3_driver_active(const A3 *a);
double a3_driver_field(const A3 *a);
double a3_driver_gradient(const A3 *a);
void a3_potential_create(A3 *a);
void a3_potential_free(A3 *a);
void a3_potential(A3 *a);
void a3_potential_seed(A3 *a, const double *phi);
void a3_transport(A3 *a, double *dn, double *dm);
int a3_geometry_cache(A3 *a, TubeGrid *g, Config3 c, int write);
int a3_geometry_identity(A3 *a, const char code[128]);
int a3_adapt(A3 *a, double **n, double **m, const unsigned char *mask, double time, int step,
             int exact_initial);
void a3_checkpoint_write(A3 *a, const double *n, const double *m, const double *marker, double time,
                         int step, double mass0, const int *quiet, double next_regrid,
                         double travel, int retries);
void a3_segment_metrics(A3 *a);
void a3_parallel_geometry(A3 *a);
void a3_face_topology(A3 *a);
void a3_checkpoint2_write(A3 *a, const double *n, const double *m, const double *marker,
                          double time, int step, double mass0, double next_regrid, double travel,
                          int retries);
int a3_checkpoint2_read(A3 *a, TubeGrid *g, Config3 c, int rank, int size, double **n, double **m,
                        double *time, int *step, double *mass0, int **quiet, double *next_regrid,
                        double *travel, int *retries);
int a3_checkpoint_read(A3 *a, TubeGrid *g, Config3 c, int rank, int size, double **n, double **m,
                       double *time, int *step, double *mass0, int **quiet, double *next_regrid,
                       double *travel, int *retries);
void a3_night_create(A3 *a, double time);
void a3_night_free(A3 *a);
void a3_night_refresh(A3 *a);
void a3_night_gather(A3 *a);
void a3_night_repartition(A3 *a, const int *cuts);
void a3_partition_plan(A3 *a, const unsigned char *work_mask, int *partition);
void a3_night_begin(A3 *a);
void a3_night_restore(A3 *a);
void a3_night_step(A3 *a, double *n, double *m, double dt);
double a3_night_loss(const A3 *a);
void a3_night_snapshot(A3 *a, const double *n, double time, int frame, double mass0);
void a3_night_electrons(A3 *a);
double a3_night_electron(A3 *a, int level, int i, int j, int k);
void a3_night_transport(A3 *a);
void a3_night_transport_begin(A3 *a);
int a3_night_predict(A3 *a, double dt);
int a3_night_correct(A3 *a, double dt);
void a3_night_accept(A3 *a);
void a3_rebase(A3 *a);
int a3_nx(const A3 *a, int level);
int a3_ny(const A3 *a, int level);
int a3_scale(int level);
int a3_macro_index(const A3 *a, int b);
int a3_macro_primary(const A3 *a, int b);
void a3_macro_init(A3 *a, A3 *reuse, int mode);
void a3_macro_prepare(A3 *a);
void a3_macro_metrics(A3 *a);
void a3_macro_tags(A3 *a, const unsigned char *fine_mask);
double a3_macro_sample(A3 *a, int field, int i, int j, int k);
void a3_macro_output(A3 *a, int frame, double time);
void a3_macro_maps(A3 *a);
#endif
