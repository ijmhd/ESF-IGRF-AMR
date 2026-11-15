/*
 * Command-line parsing and dispatch. The executable defaults are deliberately distinguished from
 * fully specified example cases; argument validation precedes model construction.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "solver.h"
#ifdef ESF_PETSC
#include <petscksp.h>
#endif
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _OPENMP
#include <omp.h>
#endif
/**
 * Parse one finite numeric argument with complete trailing-character/error checks.
 */
static double number(Solver3 *s, const char *v) {
    char *end;
    errno = 0;
    double x = strtod(v, &end);
    if (errno || end == v || *end || !isfinite(x))
        fail3(s, "invalid numeric argument");
    return x;
}
/**
 * Parse a nonnegative integer-valued argument in the supported command-line range.
 */
static int integer(Solver3 *s, const char *v) {
    double x = number(s, v);
    if (x != floor(x) || x < 0 || x > 1000000)
        fail3(s, "invalid integer argument");
    return (int)x;
}
/**
 * Create the output directory on rank zero and reject a directory that already contains run
 * metadata; never overwrite an existing run.
 */
static void directory(Solver3 *s, const char *out) {
    if (s->rank)
        return;
    char path[1200];
    snprintf(path, sizeof(path), "%s", out);
    for (char *p = path + 1; *p; p++)
        if (*p == '/') {
            *p = 0;
            if (mkdir(path, 0755) && errno != EEXIST)
                fail3(s, "mkdir failed");
            *p = '/';
        }
    if (mkdir(path, 0755) && errno != EEXIST)
        fail3(s, "mkdir failed");
    snprintf(path, sizeof(path), "%s/metadata.json", out);
    FILE *f = fopen(path, "r");
    if (f) {
        fclose(f);
        fail3(s, "output already contains a run");
    }
}
/**
 * Initialize MPI/PETSc, verify the thread support needed by OpenMP, parse/validate all CLI
 * controls, load external driver metadata and dispatch to static or AMR evolution.
 */
int main(int argc, char **argv) {
    int rank = 0, size = 1;
#ifdef ESF_MPI
#ifdef ESF_PETSC
    int petsc_argc = 1;
    char *petsc_args[] = {argv[0], NULL};
    char **petsc_argv = petsc_args;
    if (PetscInitialize(&petsc_argc, &petsc_argv, NULL, NULL))
        return 1;
#else
    MPI_Init(&argc, &argv);
#endif
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#ifdef _OPENMP
    int thread_level = 0;
    MPI_Query_thread(&thread_level);
    if (omp_get_max_threads() > 1 && thread_level < MPI_THREAD_FUNNELED) {
        if (!rank)
            fputs("OpenMP requires MPI_THREAD_FUNNELED or stronger\n", stderr);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (!rank)
        printf("Hybrid configuration: MPI=%d OpenMP=%d MPI_thread_level=%d\n", size,
               omp_get_max_threads(), thread_level);
#endif
#endif
    Solver3 s = {.rank = rank};
    Config3 c = {.night_reactions = 1,
                 .grid = {.nx = 24,
                          .ny = 40,
                          .ns = 64,
                          .span = 2.5,
                          .hmin = 250e3,
                          .hmax = 1600e3,
                          .hbase = 90e3,
                          .trace_step = 10e3},
                 .amr_levels = 2,
                 .amr_structure = 0,
                 .amr_slab = 0,
                 .amr_seed_guard = 1,
                 .amr_curve_on = .02,
                 .amr_curve_off = .01,
                 .amr_velocity_on = .1,
                 .amr_sigma_on = .05,
                 .amr_quiet_s = 60,
                 .amr_seed_amplitude = 1e-3,
                 .amr_mode = -1,
                 .amr_coo = 1,
                 .amr_work_weight = .5,
                 .amr_buffer = 2,
                 .amr_hold = 2,
                 .amr_interval = 100,
                 .amr_deplete_on = .95,
                 .amr_deplete_off = .98,
                 .amr_grad_on = .08,
                 .amr_grad_off = .04,
                 .pdm_beta = 1,
                 .degree = 13,
                 .parallel = 1,
                 .parallel_order = 1,
                 .hall = 1,
                 .regional_x = 1,
                 .max_steps = 200000,
                 .dipole_equator_nT = 31200,
                 .year = 2025.,
                 .end = 3600,
                 .output = 600,
                 .amplitude = .1,
                 .seed_height = 350e3,
                 .seed_width = 70e3,
                 .e_peak = 1e10,
                 .temperature = 2500,
                 .cfl = .35,
                 .dtmax = 2.,
                 .rtol = 1e-9};
    strcpy(c.coeff, IGRF_COEFF_FILE);
    strcpy(c.out, "igrf_run");
    int lon_set = 0, exb_set = 0;
    double seed_lon = NAN;
    for (int a = 1; a < argc; a++) {
        const char *key = argv[a];
        if (!strcmp(key, "--help")) {
            if (!rank)
                puts(
                    "esf3d --field igrf|dipole|aligned-dipole --year 2025 --nx 24 --ny 40 --ns 64\n --span-deg 2.5 --lon-deg AUTO_HK_APEX --apex-min-km 250 --apex-max-km 1600\n --base-km 90 --trace-step-km 10 --end 3600 --output-every 600\n --amplitude .1 --seed-height-km 350 --seed-width-km 70 --seed-east-km 0\n --E-peak 1e10 --temperature-sum 2500 --parallel 1 --hall 1\n --background-profile single|telf --dipole-equator-nT 31200\n --periodic-x 0 --open-apex 0 --open-ends 0 --parallel-order 1|2\n --amr-mode none|fixed|fine|dynamic --amr-case esf|translation|potential|remap\n --exb-scheme muscl|pdm4 --pdm-beta 1\n --parallel-flux rusanov|relaxation --petsc-potential 0|1\n --electric-driver FILE (CSV time_s,east_V_m; AMR, nonperiodic longitude)\n --background-east-mVm 0 --regional-x 1 (regional open longitude, also at zero imposed E)\n --cfl .35 --dt-max 2 --rtol 1e-9 --max-steps 200000 --progress-steps 0\n --out DIRECTORY --coeff FILE\n --seed-lon-deg APEX_LONGITUDE (defaults to domain center)");
#ifdef ESF_MPI
#ifdef ESF_PETSC
            PetscFinalize();
#else
            MPI_Finalize();
#endif
#endif
            return 0;
        }
        if (a + 1 >= argc)
            fail3(&s, "missing argument");
        const char *v = argv[++a];
        if (!strcmp(key, "--electric-driver")) {
            if (strlen(v) >= 1024)
                fail3(&s, "driver path too long");
            strcpy(c.driver_path, v);
            continue;
        }
        if (!strcmp(key, "--neutral-profile") || !strcmp(key, "--night-restart")) {
            char *target = !strcmp(key, "--neutral-profile") ? c.neutral_profile : c.night_restart;
            if (strlen(v) >= 1024)
                fail3(&s, "path too long");
            strcpy(target, v);
            continue;
        }
        if (!strcmp(key, "--out") || !strcmp(key, "--coeff") ||
            !strcmp(key, "--collision-profile") || !strcmp(key, "--restart") ||
            !strcmp(key, "--geometry-cache")) {
            char *target =
                !strcmp(key, "--out")
                    ? c.out
                    : (!strcmp(key, "--coeff")
                           ? c.coeff
                           : (!strcmp(key, "--restart")
                                  ? c.restart
                                  : (!strcmp(key, "--geometry-cache") ? c.geometry_cache
                                                                      : c.collision_profile)));
            if (strlen(v) >= 1024)
                fail3(&s, "path too long");
            strcpy(target, v);
            continue;
        }
        if (!strcmp(key, "--field")) {
            c.aligned_dipole = 0;
            if (!strcmp(v, "igrf"))
                c.degree = 13;
            else if (!strcmp(v, "dipole"))
                c.degree = 1;
            else if (!strcmp(v, "aligned-dipole")) {
                c.degree = 1;
                c.aligned_dipole = 1;
            } else
                fail3(&s, "unknown field");
            continue;
        }
        if (!strcmp(key, "--background-profile")) {
            if (!strcmp(v, "single"))
                c.background_profile = 0;
            else if (!strcmp(v, "telf"))
                c.background_profile = 1;
            else
                fail3(&s, "unknown initial profile");
            continue;
        }
        if (!strcmp(key, "--parallel-flux")) {
            if (!strcmp(v, "rusanov"))
                c.parallel_relaxation = 0;
            else if (!strcmp(v, "relaxation"))
                c.parallel_relaxation = 1;
            else
                fail3(&s, "unknown parallel flux");
            continue;
        }
        if (!strcmp(key, "--exb-scheme")) {
            exb_set = 1;
            if (!strcmp(v, "pdm4"))
                c.exb_pdm4 = 1;
            else if (!strcmp(v, "muscl"))
                c.exb_pdm4 = 0;
            else
                fail3(&s, "unknown ExB reconstruction");
            continue;
        }
        if (!strcmp(key, "--amr-mode")) {
            if (!strcmp(v, "none"))
                c.amr_mode = 0;
            else if (!strcmp(v, "fixed"))
                c.amr_mode = 1;
            else if (!strcmp(v, "fine"))
                c.amr_mode = 2;
            else if (!strcmp(v, "dynamic"))
                c.amr_mode = 3;
            else
                fail3(&s, "unknown AMR mode");
            continue;
        }
        if (!strcmp(key, "--amr-case")) {
            if (!strcmp(v, "esf"))
                c.amr_case = 0;
            else if (!strcmp(v, "translation"))
                c.amr_case = 1;
            else if (!strcmp(v, "potential"))
                c.amr_case = 2;
            else if (!strcmp(v, "remap"))
                c.amr_case = 3;
            else
                fail3(&s, "unknown AMR case");
            continue;
        }
        if (!strcmp(key, "--amr-criteria")) {
            if (!strcmp(v, "structure"))
                c.amr_structure = 1;
            else if (!strcmp(v, "legacy"))
                c.amr_structure = 0;
            else
                fail3(&s, "unknown AMR criterion");
            continue;
        }
        if (!strcmp(key, "--amr-density-indicator")) {
            if (!strcmp(v, "represented"))
                c.amr_residual = 1;
            else if (!strcmp(v, "total"))
                c.amr_residual = 0;
            else
                fail3(&s, "unknown density indicator");
            continue;
        }
        if (!strcmp(key, "--amr-fixture")) {
            const char *names[] = {"none", "flat", "wall", "tip", "fork", "weak"};
            c.amr_fixture = -1;
            for (int q = 0; q < 6; q++)
                if (!strcmp(v, names[q]))
                    c.amr_fixture = q;
            if (c.amr_fixture < 0)
                fail3(&s, "unknown tagging fixture");
            continue;
        }
        double x = number(&s, v);
        if (!strcmp(key, "--geometry-only"))
            c.geometry_only = integer(&s, v);
        else if (!strcmp(key, "--amr-levels"))
            c.amr_levels = integer(&s, v);
        else if (!strcmp(key, "--nx"))
            c.grid.nx = integer(&s, v);
        else if (!strcmp(key, "--night-chemistry"))
            c.night_chemistry = integer(&s, v);
        else if (!strcmp(key, "--night-reactions"))
            c.night_reactions = integer(&s, v);
        else if (!strcmp(key, "--utc-start-unix"))
            c.utc_start = x;
        else if (!strcmp(key, "--chem-Ti"))
            c.chem_Ti = x;
        else if (!strcmp(key, "--chem-Te"))
            c.chem_Te = x;
        else if (!strcmp(key, "--ny"))
            c.grid.ny = integer(&s, v);
        else if (!strcmp(key, "--ns"))
            c.grid.ns = integer(&s, v);
        else if (!strcmp(key, "--parallel"))
            c.parallel = integer(&s, v);
        else if (!strcmp(key, "--petsc-potential"))
            c.petsc_potential = integer(&s, v);
        else if (!strcmp(key, "--hall"))
            c.hall = integer(&s, v);
        else if (!strcmp(key, "--regional-x"))
            c.regional_x = integer(&s, v);
        else if (!strcmp(key, "--periodic-x"))
            c.periodic_x = integer(&s, v);
        else if (!strcmp(key, "--open-apex"))
            c.open_apex = integer(&s, v);
        else if (!strcmp(key, "--open-ends"))
            c.open_ends = integer(&s, v);
        else if (!strcmp(key, "--parallel-order"))
            c.parallel_order = integer(&s, v);
        else if (!strcmp(key, "--progress-steps"))
            c.progress_steps = integer(&s, v);
        else if (!strcmp(key, "--dipole-equator-nT"))
            c.dipole_equator_nT = x;
        else if (!strcmp(key, "--max-steps"))
            c.max_steps = integer(&s, v);
        else if (!strcmp(key, "--amr-buffer"))
            c.amr_buffer = integer(&s, v);
        else if (!strcmp(key, "--amr-slab-cells"))
            c.amr_slab = integer(&s, v);
        else if (!strcmp(key, "--amr-seed-guard"))
            c.amr_seed_guard = integer(&s, v);
        else if (!strcmp(key, "--amr-curve-on"))
            c.amr_curve_on = x;
        else if (!strcmp(key, "--amr-curve-off"))
            c.amr_curve_off = x;
        else if (!strcmp(key, "--amr-velocity-on"))
            c.amr_velocity_on = x;
        else if (!strcmp(key, "--amr-sigma-on"))
            c.amr_sigma_on = x;
        else if (!strcmp(key, "--amr-quiet-seconds"))
            c.amr_quiet_s = x;
        else if (!strcmp(key, "--amr-seed-wavelength-km"))
            c.amr_seed_lambda = x * 1000;
        else if (!strcmp(key, "--amr-branch-wavelength-km"))
            c.amr_branch_lambda = x * 1000;
        else if (!strcmp(key, "--amr-seed-amplitude"))
            c.amr_seed_amplitude = x;
        else if (!strcmp(key, "--amr-coo"))
            c.amr_coo = integer(&s, v);
        else if (!strcmp(key, "--amr-work-weight"))
            c.amr_work_weight = x;
        else if (!strcmp(key, "--amr-hold"))
            c.amr_hold = integer(&s, v);
        else if (!strcmp(key, "--amr-interval"))
            c.amr_interval = x;
        else if (!strcmp(key, "--amr-deplete-on"))
            c.amr_deplete_on = x;
        else if (!strcmp(key, "--amr-deplete-off"))
            c.amr_deplete_off = x;
        else if (!strcmp(key, "--amr-gradient-on"))
            c.amr_grad_on = x;
        else if (!strcmp(key, "--amr-gradient-off"))
            c.amr_grad_off = x;
        else if (!strcmp(key, "--pdm-beta"))
            c.pdm_beta = x;
        else if (!strcmp(key, "--year"))
            c.year = x;
        else if (!strcmp(key, "--end"))
            c.end = x;
        else if (!strcmp(key, "--output-every"))
            c.output = x;
        else if (!strcmp(key, "--span-deg"))
            c.grid.span = x;
        else if (!strcmp(key, "--lon-deg")) {
            c.grid.lon = x;
            lon_set = 1;
        } else if (!strcmp(key, "--apex-min-km"))
            c.grid.hmin = x * 1000;
        else if (!strcmp(key, "--apex-max-km"))
            c.grid.hmax = x * 1000;
        else if (!strcmp(key, "--base-km"))
            c.grid.hbase = x * 1000;
        else if (!strcmp(key, "--trace-step-km"))
            c.grid.trace_step = x * 1000;
        else if (!strcmp(key, "--amplitude"))
            c.amplitude = x;
        else if (!strcmp(key, "--seed-height-km"))
            c.seed_height = x * 1000;
        else if (!strcmp(key, "--seed-width-km"))
            c.seed_width = x * 1000;
        else if (!strcmp(key, "--seed-east-km"))
            c.seed_east = x * 1000;
        else if (!strcmp(key, "--seed-lon-deg"))
            seed_lon = x;
        else if (!strcmp(key, "--E-peak"))
            c.e_peak = x;
        else if (!strcmp(key, "--temperature-sum"))
            c.temperature = x;
        else if (!strcmp(key, "--background-east-mVm"))
            c.background_east = x * 1e-3;
        else if (!strcmp(key, "--cfl"))
            c.cfl = x;
        else if (!strcmp(key, "--dt-max"))
            c.dtmax = x;
        else if (!strcmp(key, "--rtol"))
            c.rtol = x;
        else
            fail3(&s, "unknown option");
    }
    if (c.grid.nx < 4 || c.grid.nx > 512 || c.grid.ny < 4 || c.grid.ny > 512 || c.grid.ns < 16 ||
        c.grid.ns > 1024 || (size_t)c.grid.nx * c.grid.ny * c.grid.ns > 8000000 ||
        c.grid.span <= 0 || c.grid.span > 20 || c.grid.hbase < 80e3 || c.grid.hbase > 150e3 ||
        c.grid.hmin <= c.grid.hbase || c.grid.hmax <= c.grid.hmin || c.grid.hmax > 3000e3 ||
        c.grid.trace_step < 100 || c.grid.trace_step > 20e3 || c.year < 1900 || c.year > 2030 ||
        c.end < 0 || c.output <= 0 || fabs(c.amplitude) >= 1 || c.seed_width <= 0 ||
        c.seed_east < 0 || c.e_peak < 0 || c.temperature <= 0 || c.pdm_beta < 0 ||
        c.pdm_beta > 10 || c.cfl <= 0 || c.cfl > .45 || c.dtmax <= 0 || c.rtol <= 0 ||
        c.rtol >= .01 || c.parallel > 1 || c.hall > 1 || c.regional_x > 1 || c.periodic_x > 1 ||
        c.open_apex > 1 || c.open_ends > 1 || c.parallel_order < 1 || c.parallel_order > 2 ||
        c.dipole_equator_nT <= 0)
        fail3(&s, "configuration outside supported range");
#ifndef ESF_PETSC
    if (c.petsc_potential)
        fail3(&s, "this binary has no PETSc backend");
#endif
    if (c.petsc_potential > 1)
        fail3(&s, "invalid PETSc option");
    if (c.amr_levels < 2 || c.amr_levels > 3)
        fail3(&s, "AMR supports two or three physical levels");
    if (c.amr_levels == 3 && (!c.amr_structure || !c.amr_slab))
        fail3(&s, "three levels require segmented structure AMR");
    if (c.amr_coo > 1 || c.amr_work_weight < 0 || c.amr_work_weight > 16 || c.amr_buffer > 16 ||
        c.amr_hold < 1 || c.amr_hold > 100 || c.amr_interval <= 0 || c.amr_deplete_on <= 0 ||
        c.amr_deplete_on >= c.amr_deplete_off || c.amr_deplete_off >= 1 || c.amr_grad_off <= 0 ||
        c.amr_grad_on < c.amr_grad_off)
        fail3(&s, "invalid AMR refinement controls");
    if (c.periodic_x && (!c.aligned_dipole || c.background_east != 0 || c.seed_east > 0))
        fail3(
            &s,
            "periodic longitude currently requires an aligned dipole, zero background E and a periodic seed");
    if (!lon_set) {
        Igrf m;
        V3 apex;
        if (solver3_load_field(&c, &m) ||
            tube_apex(&m, geodetic_position(22.3, 114.2, 350e3), &apex))
            fail3(&s, "cannot locate target flux tube");
        c.grid.lon = atan2(apex.y, apex.x) * 180 / 3.14159265358979323846;
    }
    if (isfinite(seed_lon)) {
        c.seed_offset_deg = seed_lon - c.grid.lon;
        if (fabs(c.seed_offset_deg) >= .5 * c.grid.span)
            fail3(&s, "seed apex longitude must lie strictly inside the domain");
        if (c.periodic_x)
            fail3(&s, "an explicit seed longitude requires a regional domain");
    }
    if (c.amr_slab && (!c.amr_structure || c.amr_slab > c.grid.ns || c.grid.ns % c.amr_slab))
        fail3(&s, "segmented density requires structure criterion and ns divisible by slab cells");
    if (c.amr_seed_guard > 1 || c.amr_curve_off <= 0 || c.amr_curve_on <= c.amr_curve_off ||
        c.amr_velocity_on <= 0 || c.amr_sigma_on <= 0 || c.amr_quiet_s < 0 ||
        c.amr_seed_amplitude <= 0 || c.amr_seed_amplitude >= 1 || c.amr_seed_lambda < 0 ||
        c.amr_branch_lambda < 0)
        fail3(&s, "invalid structure refinement controls");
    ElectricDriver *driver = NULL;
    if (c.driver_path[0]) {
        if (c.amr_mode < 0 || c.background_east != 0)
            fail3(&s, "CSV driver requires AMR backend and excludes constant background-east");
        driver = calloc(1, sizeof(*driver));
        if (!driver || electric_driver_load(driver, c.driver_path, c.end))
            fail3(&s,
                  "invalid electric driver: require ordered time_s,east_V_m with coverage0..end");
        c.driver = driver;
    }
    if (c.amr_mode >= 0) {
        if (c.neutral_profile[0]) {
            FILE *nf = fopen(c.neutral_profile, "rb");
            char magic[8];
            if (!nf)
                fail3(&s, "neutral profile open");
            c.spatial_neutral = fread(magic, 1, 8, nf) == 8 && !memcmp(magic, "ESFNEU01", 8);
            fclose(nf);
            if (c.spatial_neutral)
                c.neutral_grid_hash = solver3_file_hash(&s, c.neutral_profile);
        }
        if (c.night_chemistry) {
            if (c.night_chemistry != 1 || c.night_reactions > 1 || !c.parallel_relaxation ||
                !c.amr_structure || !c.amr_slab || c.amr_case || !c.neutral_profile[0] ||
                c.utc_start < 946684800 || c.utc_start > 4102444800 || c.chem_Ti < 200 ||
                c.chem_Te < 200 || fabs(c.chem_Ti + c.chem_Te - c.temperature) > 1e-10 ||
                !c.e_peak || (c.restart[0] && !c.night_restart[0]))
                fail3(
                    &s,
                    "night closure requires physical segmented AMR, neutral table, UTC, consistent Ti+Te, positive initial E layer and paired night restart");
        } else if (c.night_restart[0] || c.neutral_profile[0] || c.utc_start || c.chem_Ti ||
                   c.chem_Te)
            fail3(&s, "night inputs require --night-chemistry1");
        if (!exb_set)
            c.exb_pdm4 = 1;
#ifndef ESF_PETSC
        fail3(&s, "the composite AMR potential requires an ESF_PETSC build");
#else
        directory(&s, c.out);
        int result = amr3_run(c, rank, size);
        free(driver);
        PetscFinalize();
        return result;
#endif
    }
    if (c.night_chemistry)
        fail3(&s, "night chemistry is supported only on segmented AMR, not the static backend");
    directory(&s, c.out);
    solver3_create(&s, c, rank, size);
    double time = 0;
    int step = 0;
    if (c.restart[0])
        solver3_restart(&s, c.restart, &time, &step);
    double start_time = time;
    int start_step = step;
    if (!rank && c.restart[0]) {
        char path[1200];
        snprintf(path, sizeof(path), "%s/restart_info.json", c.out);
        FILE *f = fopen(path, "w");
        if (!f)
            fail3(&s, "cannot record restart");
        fprintf(
            f,
            "{\"time_s\":%.17g,\"step\":%d,\"original_mass\":%.17g,\"integrated_boundary_particle_flux\":%.17g,\"potential_warm_start\":true}\n",
            time, step, s.mass0, s.boundary_total);
        fclose(f);
    }
    double evolution_start = 0, evolution_wall = 0, timing_end[PROF_COUNT] = {0};
    long iterations_end = 0;
    double next = time;
    int frame = 0, start = s.first * c.grid.ns, end = (s.first + s.count) * c.grid.ns;
    for (;;) {
        if (time >= next - 1e-9 || time >= c.end - 1e-9) {
            if (frame > 0 && time >= c.end - 1e-9) {
                evolution_wall = wall3() - evolution_start;
                memcpy(timing_end, s.timing, sizeof(s.timing));
                iterations_end = s.potential_iterations;
            }
            if (solver3_rhs(&s, s.rho, s.momentum, s.rhs_n, s.rhs_m))
                fail3(&s, "potential solve failed at output");
            solver3_output(&s, time, frame++, step);
            if (frame == 1) {
#ifdef ESF_MPI
                MPI_Barrier(MPI_COMM_WORLD);
#endif
                memset(s.timing, 0, sizeof(s.timing));
                s.potential_iterations = 0;
                evolution_start = wall3();
            } else if (time >= c.end - 1e-9) {
                memcpy(s.timing, timing_end, sizeof(s.timing));
                s.potential_iterations = iterations_end;
            }
            next += c.output;
            if (time >= c.end - 1e-9)
                break;
        }
        if (step >= c.max_steps)
            fail3(&s, "step limit");
        /* Leave room for a changing accepted-state wave speed. Without this
         * margin, a tiny stage-speed increase repeatedly halves CFL-limited
         * steps and doubles potential solves in fine-grid studies. */
        double dt = fmin(c.dtmax, .9 * c.cfl / fmax(s.cfl_rate, 1e-30));
        dt = fmin(dt, fmin(c.end - time, next - time));
        memcpy(s.old_n + start, s.rho + start, s.n * sizeof(double));
        memcpy(s.old_m + start, s.momentum + start, s.n * sizeof(double));
        int ok = 0;
        double boundary_increment = 0;
        for (int attempt = 0; attempt < 24; attempt++) {
            for (int k = start; k < end; k++) {
                s.rho[k] = s.old_n[k];
                s.momentum[k] = s.old_m[k];
            }
            if (solver3_rhs(&s, s.rho, s.momentum, s.rhs_n, s.rhs_m))
                fail3(&s, "stage potential did not converge");
            dt = fmin(dt, .9 * c.cfl / fmax(s.cfl_rate, 1e-30));
            double boundary0 = s.boundary_rate;
            if (dt * s.cfl_rate <= c.cfl * (1 + 1e-12)) {
                for (int k = start; k < end; k++) {
                    double decay, p1, p2;
                    exponential_weights(s.g.cell[k].nu, dt, &decay, &p1, &p2);
                    s.first_rhs_m[k] = s.rhs_m[k];
                    s.stage_n[k] = s.rho[k] + dt * s.rhs_n[k];
                    s.stage_m[k] = decay * s.old_m[k] + dt * p1 * s.rhs_m[k];
                }
                if (solver3_valid(&s, s.stage_n, s.stage_m)) {
                    if (solver3_rhs(&s, s.stage_n, s.stage_m, s.rhs_n, s.rhs_m))
                        fail3(&s, "second stage potential did not converge");
                    double boundary1 = s.boundary_rate;
                    if (dt * s.cfl_rate <= c.cfl * (1 + 1e-12)) {
                        for (int k = start; k < end; k++) {
                            double decay, p1, p2;
                            exponential_weights(s.g.cell[k].nu, dt, &decay, &p1, &p2);
                            s.rho[k] = .5 * s.old_n[k] + .5 * (s.stage_n[k] + dt * s.rhs_n[k]);
                            s.momentum[k] = decay * s.old_m[k] +
                                            dt * ((p1 - p2) * s.first_rhs_m[k] + p2 * s.rhs_m[k]);
                        }
                        if (solver3_valid(&s, s.rho, s.momentum)) {
                            boundary_increment = .5 * dt * (boundary0 + boundary1);
                            ok = 1;
                            break;
                        }
                    }
                }
            }
            dt *= .5;
            s.retries++;
        }
        if (!ok || dt < 1e-12)
            fail3(&s, "positive/CFL-compliant step not found");
        s.boundary_total += boundary_increment;
        time += dt;
        step++;
        if (!rank && c.progress_steps && step % c.progress_steps == 0) {
            fprintf(stderr, "progress step=%d time_s=%.9g dt=%.3g retries=%d wall_s=%.2f\n", step,
                    time, dt, s.retries, wall3() - s.wall0);
            fflush(stderr);
        }
    }
    if (!rank) {
        char path[1200];
        snprintf(path, sizeof(path), "%s/completion.json", c.out);
        FILE *f = fopen(path, "w");
        if (!f)
            fail3(&s, "cannot write completion");
        fprintf(
            f,
            "{\"time_s\":%.17g,\"start_time_s\":%.17g,\"steps\":%d,\"segment_steps\":%d,\"solves\":%d,\"retries\":%d,\"maximum_potential_relative_residual\":%.17g,\"wall_s\":%.6f}\n",
            time, start_time, step, step - start_step, s.solves, s.retries,
            s.max_potential_residual, wall3() - s.wall0);
        fclose(f);
    }
#ifdef ESF_MPI
    MPI_Allreduce(MPI_IN_PLACE, &evolution_wall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    double *all = calloc((size_t)size * PROF_COUNT, sizeof(double));
#ifdef ESF_MPI
    MPI_Gather(s.timing, PROF_COUNT, MPI_DOUBLE, all, PROF_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    memcpy(all, s.timing, sizeof(s.timing));
#endif
    if (!rank) {
        const char *names[] = {
            "state_exchange",     "column_integration", "coefficient_exchange",
            "potential_assembly", "potential_solve",    "potential_broadcast_and_wait",
            "face_flux",          "cell_rhs",           "cfl_and_boundary_reductions",
            "state_validation"};
        char path[1200];
        snprintf(path, sizeof(path), "%s/profile.json", c.out);
        FILE *f = fopen(path, "w");
        if (!f)
            fail3(&s, "profile output failed");
        fprintf(
            f,
            "{\"steps\":%d,\"compute_only_wall_s\":%.9g,\"potential_iterations\":%ld,\"ranks\":[",
            step, evolution_wall, s.potential_iterations);
        for (int r = 0; r < size; r++) {
            fprintf(f, "%s{", r ? "," : "");
            for (int k = 0; k < PROF_COUNT; k++)
                fprintf(f, "%s\"%s\":%.9g", k ? "," : "", names[k], all[r * PROF_COUNT + k]);
            fprintf(f, "}");
        }
        fprintf(f, "]}\n");
        fclose(f);
    }
    free(all);
    solver3_free(&s);
#ifdef ESF_MPI
#ifdef ESF_PETSC
    PetscFinalize();
#else
    MPI_Finalize();
#endif
#endif
    return 0;
}
