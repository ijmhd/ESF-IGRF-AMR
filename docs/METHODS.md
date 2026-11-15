# Physics, numerical methods and call sequence

The complete derivations, constants, algorithm branches and implementation
reference are in [methodology/main.tex](../methodology/main.tex). This guide
provides a reading map tied to actual code.

## State and equations

Stored density is N=n/(1e12 m^-3), parallel momentum M=N*u_parallel. All positions,
volumes and velocities otherwise use SI. With night mode, quasineutral density
is Ne=NOplus+NNOplus+NO2plus; without it the prescribed physical molecular E
population supplies the additional electrons/conductivity.

Each transported species has conservative continuity. Its along-field momentum
uses inertia/advection, gravity, fixed-temperature ion pressure, common electron
pressure and linear neutral drag. In the relaxation branch, pressure/gravity
enter the cell source and a fitted mass flux; momentum advects with that same
mass flux and upwind velocity. This is a reduced discretization, not a full
multi-fluid energy/momentum system.

O+ perpendicular mass motion uses E-cross-B. Molecular perpendicular motion
includes species-specific collisional Hall/Pedersen electric mobility and gravity
drift. The polarization forcing also includes gravity currents. O+ mass transport
therefore uses the standard reduced F-region drift approximation rather than
the full perpendicular species force balance. Molecular perpendicular pressure
diffusion and ion-ion momentum exchange are omitted.

The electrostatic equation enforces field-line-integrated current continuity.
Its weak form retains the full nonorthogonal Pedersen/Hall tensor and gravity
current forcing. Each traced line is equipotential. Bilinear constrained basis
functions handle hanging nodes; PETSc COO tables cache only geometry/basis
algebra, never old plasma coefficients. FGMRES/GAMG solves the current system,
and an explicit infinity-norm true-residual check enforces the requested tolerance.

Boundary conditions: polarization potential is fixed at transverse apex-height
extremes; regional longitude uses `d_xi(phi_pol)=0` with Hall/geometric current
terms retained. This is not zero normal total current. Imposed eastward E is
mapped through a linear potential term at actual stage time. Open particle
boundaries use the code's extrapolated exterior state; closed along-field ends
reflect parallel velocity and have zero particle flux.

## Spatial methods

PDM4 is a four-cell physical-volume cubic reconstruction for O+ transverse
transport. Moment weights reproduce degrees0-3 on the mapped hexahedra. The
limiter clips to donor/downstream bounds, suppresses extrema and limits the
correction to `2*beta*|donor-upstream|`. beta=1 is the release case. PDM4 is not
PDM7 and does not make the complete coupled method fourth order in time.
Boundary stencils use the explicit lower-order fallback. Parallel reconstruction
and molecular donor-cell perpendicular fluxes are separate schemes.

For a parallel face, `theta=1/[1+(nu*ds/cs)^2]` blends an inertial Rusanov particle
flux with a fitted collisional flux. The Bernoulli function is evaluated stably
near zero/large arguments. Common Ne appears in the fitted electron potential
difference and cell pressure source. Exact zero density uses a bounded secant;
positive densities are not replaced by a floor.

## Complete accepted step

1. Check tag clock/travel; when masks change, migrate, conservatively remap,
   rebuild constrained field layout and solve the remapped potential.
2. Propose dt from the previous transport rate, maximum dt and stop/output/
   driver/regrid clocks. Back up O+, marker, molecules and chemical budgets.
3. Apply chemistry for dt/2. Each local reaction preserves positivity; O+ is
   volume-averaged back to its actual AMR volume.
4. `rhs`: O+ halo -> prepare restrictions/virtual states -> molecular halo and
   conductivity refresh -> common Ne -> integrated potential solve -> potential
   face derivatives -> O+ flux/source/rates -> molecular flux/source/rates.
5. If this actual post-chemistry state requires a smaller initial CFL dt, restore
   every backed-up field/budget and retry. This is counted separately from
   predictor/corrector failure retries.
6. Predict all densities and momenta together. Density uses the explicit first
   stage; momentum uses exact fixed-nu drag weights. Reject invalid states.
7. Run the entire `rhs` chain on the predicted O+ AND molecular state at t+dt.
8. Correct density with the two-stage average and momentum with ETD2. Check the
   predicted-stage CFL and all densities/momenta/marker. On failure restore and
   reduce dt; do not reuse a different state's phi or flux.
9. Commit the accepted molecular boundary increment, apply the final chemistry
   half-step, and update global boundary/travel histories and physical time.

No level subcycling is used. ETD weights may be reused only between predictor
and corrector of the exact same attempted dt and mesh. Log-pressure caches are
refreshed at every molecular RHS stage. MPI/PETSc calls occur on the main thread.

The owned Ne averaging cache has a separate validity rule: it is cleared every
time `a3_night_electrons` forms the electron state. A first pressure/parallel-flux
request performs the original volume-weighted sum and division; repeated
requests at the same level/cell reuse that result. This includes the original
finest-level `Ne*V/V` arithmetic. Outside ownership, the original direct path is
retained. Lazy writes occur on the main thread. It is not a previous-step Ne
approximation or an interpolation of an old pressure field.

## Chemistry and prescribed inputs

The four reactions are O+ + N2 -> NO+ + N, O+ + O2 -> O2+ + O, NO+ + e -> neutrals,
and O2+ + e -> neutrals. Thermal zero-drift fits are recorded from the pinned
[GEMINI3D rate source](https://github.com/gemini3d/gemini3d/blob/d70446a712d47e7617b25c033bfcace64207590f/src/sources/sources.f90).
The low-effective-temperature conversion branches reject temperatures outside
their implemented range. Symmetric exact scalar suboperators retain a charge
budget; this is not full photochemistry. Existing species keep their velocity
during a local reaction; newly populated zero-density species start with zero
parallel momentum. The fixed neutral bath supplies/removes reaction momentum.

Collision expressions follow the separately audited
[GLOW conductance source](https://github.com/gemini3d/glow/blob/master/conduct.f90).
The implemented electron collision expression has an explicit factor4. The
production O+ collision input also inherits its prescribed generator policy;
these empirical choices are inputs/approximations, not fitted away during
performance work. Exact formulas and SI/cm conversion appear in the manuscript.
The idealized built-in collision background and synthetic demo neutral table are
not an IRI/MSIS event background.

The night guard tests geometric Earth shadow of every modeled E-region cell at
1800s intervals and the requested end. No photoionization is supplied. This is
a coarse applicability guard, not precise terminator integration or an MLT model.
Temperatures, neutrals, B and geometry are prescribed, not dynamically evolved.

## Source reading map

`main.c` -> `amr3_run` -> `rhs` -> `a3_prepare` / `a3_potential` /
`a3_transport` / `a3_night_transport` is the central chain. Read
`amr3_structure.c`, `amr3_macro.c`, `amr3_adapt.c` for AMR; `pdm.c` for
reconstruction; `night.c` / `night_transport.c` for local physics kernels;
`geometry.c` / `igrf.c` for magnetic geometry; and the checkpoint modules for
portable persistence. [FUNCTIONS.md](FUNCTIONS.md) covers every runtime function.
