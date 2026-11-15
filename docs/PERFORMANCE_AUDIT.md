# Performance audit and remaining work

The completed numbers and exact qualification scope are in
[VALIDATION.md](VALIDATION.md). This audit explains which implementation costs
remain and what would have to change to remove them. It does not infer an
exclusive time budget from overlapping timers.

## Changes retained in v0.1.0

* Build three-level raw metrics once, then apply actual owner masks. A restart
  still rebuilds/remasks the required geometry; it does not reuse a stale mask.
* Reuse immutable macro volume/background tables with separate ownership.
* Reuse identical local conversion coefficients within the same reaction call,
  while reevaluating nonlinear recombination at the updated density.
* Reuse ETD coefficients between predictor/corrector of one attempted dt.
* Cache owned volume-weighted Ne samples only within a rebuilt electron state.
  The first sample retains the original floating-point sum/division.
* Give independent cell work and disjoint molecular lines to OpenMP workers.
  Keep the field-coordinate index contiguous; preserve face order on a line.
* Store real MPI status arrays and use explicit C11 array qualifications.
  Compiler warnings are fixed at their cause, not disabled.
* Expose initial CFL adjustments independently from failure retries.

The prior qualified reference already replaced repeated complete molecular
gathers at regridding with owned-intersection migration, removed unused O+
logarithms and cached current-RHS molecular pressure logs. Those gains are
included in the reference and are not counted a second time.

The final combined release is compared with that reference. Measurements of
earlier candidates use separate allocations and cannot isolate the gain from
each individual cache or loop edit. A token-parity proof applies to the initial
annotation pass; later intentional numerical-source edits are instead bound to
their own exact identity and saved-state qualification.

## Row ownership limits rank balance

The three-level science grid has only80 physical root rows. With64 MPI ranks,
48 ranks own one row and16 own two at the audited13500s state. Molecular
owned-cell maximum/mean is1.6; actual O+ activity maximum/mean is about1.658.
Adding ranks does not create finer ownership granularity. The root-row prefix
cost estimate includes canonical work, but does not exactly predict chemistry,
macro coarsening, reconstruction or ghost work.

The O+ halo extends six internal rows on each side in the three-level path.
At the audited layout the arithmetic mean of per-rank ghost/backing-owned
ratios is5.219. That is a storage/work count, not a percentage of runtime.
It explains why many very thin MPI slabs can be expensive even when their
owned cells are well balanced by a scalar estimate.

A two-dimensional tile/block partition could improve granularity and reduce
halo surface/volume. It would require a new ownership/neighbor representation,
multi-neighbor halos, arbitrary old/new patch intersections, coarse/fine face
ownership, potential-data assembly and portable-restart redistribution. The
current distributed PETSc solve is not by itself such a density decomposition.
This change requires new conservation, migration and developed-state checks;
it is not an unreviewed pragma or launch-flag change.

## O+ activity is only part of the work

At the same audited state there are2219120 active O+ volumes,4739264 backing
O+ slots and5242880 canonical cells for **each** molecular species. Geometry
and several coefficient/state arrays retain global canonical backing on every
rank. Only owned/refreshed ghost values are valid during evolution.

Consequences:

1. Fewer O+ active volumes do not reduce molecular transport/chemistry cells
   by the same factor. Refining O+ to L2 can increase canonical molecular cost
   across the whole domain even if O+ L2 occupies a small region.
2. A reconstructed display array has canonical resolution without that many
   independently evolved O+ values. Counting image pixels overstates activity.
3. Replicated geometry/backing limits memory scaling. A roughly4.7-4.9GiB
   per-process peak in the audited large configuration is not all owned plasma.
   A sum of process high-water marks is not simultaneous job RSS.
4. Blindly skipping zero-owner slots can break virtual children and PDM
   sampling. An active-index optimization must classify physical updates,
   restriction/reconstruction dependencies and output backing separately.

Distributing immutable geometry and making molecular state genuinely local
would reduce replicated memory and some traffic. Molecular AMR is a further
algorithm change: chemical charge ledgers, common Ne/electron pressure,
integrated conductivity and coupled species remaps must remain consistent
across every differing resolution. It is not already present in this release.

## Coupled attempts contain substantial repeated work

The new counters distinguish accepted steps, attempts and initial CFL
adjustments. In the qualified interval400 steps require799 attempts because
399 post-chemistry initial CFL proposals are reduced. That produces1200 RHS
calls while the legacy failure-retry count remains0.

An RHS includes state preparation, current Ne/conductivity, a potential solve
and all species fluxes. The first half-chemistry operation changes the state,
so a rejected larger-dt attempt cannot simply donate its potential or flux to
the smaller-dt attempt. Doing so would alter the coupled method. A future
proposal predictor could reduce repeated attempts, but it must respect the
post-chemistry stability condition and must be qualified as an explicit
time-integration change. This version retains the existing dt policy.

## OpenMP has a deliberate serial fraction

Molecular line and independent cell kernels have disjoint writes. O+ AMR
face traversal, lazy child reconstruction and lazy Ne sampling remain on the
main thread. MPI/PETSc also use that thread. More workers cannot accelerate
those sections directly. Molecular axis sweeps retain barriers because each
axis contributes to the same cell arrays.

Threading arbitrary O+ faces would race both on neighboring cells and on
lazy reconstruction validity/data. A safe future implementation could separate
read-only preparation from flux evaluation and use owned face buffers plus a
deterministic cell gather. Its additional buffer cost, cache behavior and
coarse/fine conservation need measurement; there is no demonstrated gain yet.

## Measurement rules

Use the complete immediate process wait for elapsed-time comparisons. Bind
the exact checkpoint pair, physics, grid, dt/solver policy, output/checkpoint
cadence, executable and worker affinity. Report active threads separately
from reserved CPUs. The qualified32x2 layout is slower than64x1, so a hybrid
launch is not automatically a recommendation.

Rank0 phase timers overlap: regrid includes preparation/solving, and the
legacy chemistry timer includes response refresh and halo wait. Separate
maxima cannot be subtracted or summed into100%. PMPI instrumentation can
explain collectives but changes the observation setup; an observed reference
must be paired with an observed candidate. A short dense-output run cannot
predict a20000s production wall time or establish asymptotic strong scaling.

Further optimization is justified, particularly ownership and canonical
backing. It is a substantive next version with its own numerical evidence,
not a claimed completed feature of the frozen student release.
