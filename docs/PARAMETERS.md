# Complete executable parameter reference

Defaults below are from `src/main.c`, not the launcher or production case. All
options except `--help` require one value. Units shown here are CLI units;
lengths become metres internally. Read the interaction constraints in each row.

| Option | Type / units | Executable default | Meaning / interactions |
|---|---|---|---|
| `--help` | flag | n/a | Print the short help and exit; this full reference includes extension options absent from short help. |
| `--field` | enum | igrf | igrf: degree13; dipole: actual degree1 Gauss coefficients; aligned-dipole: axis-aligned fixture with specified equatorial strength. |
| `--year` | decimal year | 2025 | Coefficient interpolation epoch, 1900-2030; independent of the chemistry UTC clock. |
| `--coeff` | path | build-source data/igrf14coeffs.txt | Official IGRF-14 coefficient file; pass a repository-relative or absolute path explicitly for relocatable runs. |
| `--lon-deg` | degrees east | trace HK350km apex | Geographic APEX longitude at domain centre, not geographic longitude at HK350km. Automatic target is geodetic22.3N114.2E350km. |
| `--span-deg` | degrees | 2.5 | Total apex-longitude width, positive and <=20. |
| `--nx` | integer | 24 | User physical L0 longitude cells; >=4 and <=512, subject to canonical size limits. |
| `--ny` | integer | 40 | User physical L0 apex-height cells; >=4 and <=512. Three-level MPI rank count cannot exceed this value. |
| `--ns` | integer | 64 | Whole-conjugate along-field segments, 16-1024; never refined by transverse AMR. |
| `--apex-min-km` | km | 250 | Lowest tube apex spherical altitude; must exceed base altitude. |
| `--apex-max-km` | km | 1600 | Highest tube apex spherical altitude, above minimum and <=3000. |
| `--base-km` | km | 90 | Both conjugate endpoints at spherical altitude above6371.2km radius; 80-150. |
| `--trace-step-km` | km | 10 | Maximum RK4 tracing step, 0.1-20; tracing also reduces steps near the base. |
| `--geometry-only` | 0/1 | 0 | Construct/load geometry and return before plasma evolution on AMR path. |
| `--geometry-cache` | path | none | Native-ABI canonical geometry cache; exact parameter/input/code/ABI matching required. |
| `--background-profile` | enum | single | single: Chapman peak438km scale70km; telf: sum of O+ Chapman populations with peaks438/600km, scales80/200km, relative amplitudes1/.01. |
| `--E-peak` | m^-3 | 1e10 | Initial physical molecular E population, peak110km scale12km. Fixed without night module; initially NO+ with night module. |
| `--amplitude` | fraction | 0.1 | Seed fractional amplitude; absolute value <1. Positive values make a depletion in the configured envelope. |
| `--seed-height-km` | km | 350 | Spherical-altitude seed envelope centre. |
| `--seed-width-km` | km | 70 | Positive vertical seed envelope width; see exact initialization equation in methods. |
| `--seed-east-km` | km | 0 | Nonnegative horizontal seed width; zero selects the alternate wavelength pattern, not an eastward displacement. |
| `--seed-lon-deg` | apex degrees | domain centre | Seed apex longitude, strictly inside the regional domain; explicit setting is incompatible with periodic-x. |
| `--temperature-sum` | K | 2500 | Fixed Ti+Te used by O+ parallel pressure/sound speed. Night module requires exact consistency with chem-Ti+chem-Te. |
| `--dipole-equator-nT` | nT | 31200 | Aligned-dipole equatorial amplitude; positive; does not replace all coefficients in ordinary degree-one dipole mode. |
| `--parallel` | 0/1 | 1 | Include field-aligned transport/momentum. Zero is a deliberately reduced physical case. |
| `--parallel-order` | 1 or2 | 1 | Order of the bounded along-field face reconstruction; distinct from transverse PDM4. |
| `--parallel-flux` | enum | rusanov | rusanov or relaxation. Night module requires relaxation; includes collision/pressure/gravity fitted mass flux. |
| `--exb-scheme` | enum | muscl static; pdm4 AMR | muscl or pdm4; explicit selection is recommended. PDM4 is a four-cell volume-moment cubic stencil with limiter. |
| `--pdm-beta` | dimensionless | 1 | PDM limiter correction bound2*beta*upstream difference; supported0-10. beta=0 reduces to donor-cell at the limiter. |
| `--hall` | 0/1 | 1 | Retain integrated Hall conductivity cross terms; disabling changes physics. |
| `--periodic-x` | 0/1 | 0 | Longitude periodicity; physical branch requires aligned dipole, zero imposed E and no localized horizontal seed. |
| `--regional-x` | 0/1 | 1 | Regional open longitude particle boundary and polarization d_xi(phi)=0 side condition; this is NOT zero total normal current. |
| `--open-apex` | 0/1 | 0 | Open particle boundaries at transverse apex-height extremes; required for external driver cases. |
| `--open-ends` | 0/1 | 0 | Open particle boundaries at both along-field endpoints; otherwise reflective/closed. |
| `--background-east-mVm` | mV/m | 0 | Prescribed constant eastward electric field, added through mapped potential; excludes CSV driver. AMR driver requires regional nonperiodic and open-apex geometry. |
| `--electric-driver` | CSV path | none | time_s,east_V_m strictly increasing from0 with coverage through absolute end; stage interpolation and step landing on knots. |
| `--petsc-potential` | 0/1 | 0 | PETSc backend for uniform reference path. AMR inherently requires PETSc even though this static-path selector defaults0. |
| `--rtol` | dimensionless | 1e-9 | Requested true potential infinity-residual scale; positive and <.01; AMR accepts <=1.05*rtol after explicit residual computation. |
| `--cfl` | dimensionless | 0.35 | Positive and <=.45; initial proposal uses0.9*cfl/rate and stage acceptance uses cfl. Teaching/qualification cases explicitly use.25. |
| `--dt-max` | s | 2 | Positive maximum attempted step; all active levels/species share the accepted dt, with no subcycling. |
| `--end` | s | 3600 | Absolute elapsed-model stop time. On restart this is not additional duration. |
| `--output-every` | s | 600 | Positive frame/checkpoint interval; initial and final frames also written. Not a physics timescale or an AMR interval. |
| `--max-steps` | integer | 200000 | Absolute step counter limit, including history restored by restart. |
| `--progress-steps` | integer | 0 | Progress-print step interval; zero disables periodic progress lines but not essential diagnostics. |
| `--out` | directory | igrf_run | New run directory; an existing metadata.json is rejected to preserve prior outputs. |
| `--restart` | path | none | Native AMR checkpoint or explicitly supported uniform-state import; physical/numerical signatures must match. |
| `--night-restart` | path | none | Molecular ESFNIT02 paired with exact native checkpoint; required for night continuation. |
| `--amr-mode` | enum | unset = static reference | none,fixed,fine,dynamic select the AMR backend. There is no literal static enum; omitting the option selects the uniform reference path. |
| `--amr-levels` | 2 or3 | 2 | Number of PHYSICAL levels. Three requires structure criterion and nonzero along-field slab size. |
| `--amr-case` | enum | esf | esf,translation,potential,remap. Latter three are verification fixtures, not physical ESF forecasts. |
| `--amr-criteria` | enum | legacy | legacy percentile/depletion or structure local solution-based detector. Released examples explicitly select structure. |
| `--amr-density-indicator` | enum | total | total or represented; represented removes a fitted local background-shaped component for detection only, never from evolved density. |
| `--amr-slab-cells` | integer | 0 | Along-field segment grouping for density tags. Zero means column-wide legacy layout; nonzero must divide ns and requires structure. Qualification uses4. |
| `--amr-interval` | s | 100 | Maximum tag-check interval. Accumulated transverse travel can trigger earlier checks; not a fixed actual regrid period. |
| `--amr-buffer` | integer | 2 | Transverse root-cell buffer around requested structure; 0-16. Qualification uses1. |
| `--amr-work-weight` | dimensionless | 0.5 | Additive canonical-field work estimate in root-row MPI partition planner; does not change physical equations. |
| `--amr-coo` | 0/1 | 1 | Use cached constrained COO assembly; zero selects the reference element assembly algebra. |
| `--amr-gradient-on` | dimensionless | 0.08 | Refine threshold for normalized detector gradient; must be >=positive off threshold. |
| `--amr-gradient-off` | dimensionless | 0.04 | Lower keep-refined gradient threshold, providing spatial hysteresis. |
| `--amr-curve-on` | dimensionless | 0.02 | Refine threshold for normalized density Hessian norm. |
| `--amr-curve-off` | dimensionless | 0.01 | Lower keep-refined curvature threshold. |
| `--amr-velocity-on` | dimensionless | 0.1 | Threshold for interval-scaled drift curvature divided by local spacing; keep-refined threshold is half. |
| `--amr-sigma-on` | dimensionless | 0.05 | Integrated conductivity-tensor curvature threshold for FIELD requests, half when already refined. |
| `--amr-quiet-seconds` | s | 60 | Physical quiet and minimum residence time before coarsening; qualified production uses60. Zero-hold simultaneous coarsening is not qualified; see known limitations. |
| `--amr-seed-guard` | 0/1 | 1 | Advect a material seed-protection marker and prevent under-resolved protected wavelengths. Fresh perturbed structure runs need explicit positive seed wavelength when enabled. |
| `--amr-seed-wavelength-km` | km | 0 = unspecified | Shortest declared protected seed wavelength, not inferred from grid. Sampling target32 cells/wavelength; unresolved locations are reported. |
| `--amr-branch-wavelength-km` | km | 0 = disabled | Additional protected wavelength when resolved walls/branches are tagged; zero disables this requirement. |
| `--amr-seed-amplitude` | fraction | 0.001 | Declared weak-seed amplitude used in protection/information-loss tolerance; distinct from density seeding amplitude. |
| `--amr-deplete-on` | ratio | 0.95 | Legacy detector depletion refinement threshold; below off threshold and positive. |
| `--amr-deplete-off` | ratio | 0.98 | Legacy keep-refined depletion threshold, <1. Not used as a density floor. |
| `--amr-hold` | check count | 2 | Legacy quiet-check count before coarsening; structure detector instead uses quiet-seconds. |
| `--amr-fixture` | enum | none | none,flat,wall,tip,fork,weak synthetic detector fixtures; explicitly separate these from scientific seeding. |
| `--night-chemistry` | 0/1 | 0 | Enable the full transported molecular module, not just local chemistry. Requires physical segmented structure AMR, neutral input, UTC, positive E peak and consistent temperatures. |
| `--night-reactions` | 0/1 | 1 | With night module on, zero disables only reactions; NO+/O2+ still transport. It is not a fixed-reservoir switch. |
| `--neutral-profile` | path | none | Five-column height table or strict spatial ESFNEU01 grid; physical values must cover the entire canonical conjugate geometry. |
| `--collision-profile` | CSV path | none | height_m,nu_O_s1,nu_e_s1; strictly increasing, positive and fully covering geometry. Without it use documented idealized analytic rates. |
| `--utc-start-unix` | Unix UTC s | 0 = unspecified | Absolute UTC at model time0, required for night module; 2000-2100 range checked. Not HKT, LST or MLT. |
| `--chem-Ti` | K | 0 = unspecified | Fixed ion temperature for all night species, 200-5000 and applicable reaction-fit range. |
| `--chem-Te` | K | 0 = unspecified | Fixed electron temperature for night closure, 200-5000. Ti+Te must equal temperature-sum. |

Dimension limits also include `nx*ny*ns <= 8000000` before AMR processing.
On AMR, the internally expanded root grid must satisfy `nx_internal*ny_internal*ns <= 2000000`.
Three-level internal nx/ny are twice user dimensions, and canonical dimensions are twice internal dimensions.
No run parameter may be inferred from a diagnostic clock. Use the checked-in explicit examples as a starting point.

The historical defaults are retained for reproducibility; they do not imply that an otherwise unspecified combination is a qualified science case. In particular, the structure examples explicitly choose PDM4, relaxation flux, slab4, seed-guard0, CFL.25 and rtol1e-10.
