# Exact-work performance update

This development branch preserves stable tag `v0.1.0`. It changes execution
work and private caches; it does not change physical inputs, PDM4 beta1,
CFL0.25, dt maximum2s, true-potential-residual request1e-10, AMR criteria,
quiet time, accepted-step policy or output/checkpoint cadence.

The current qualified numerical source is
`03a3c037cb75d2fbc640507d5b92a584158e31160124bd6753119b0834ceee76`.
[qualification.json](../validation/qualification.json) binds its executable,
flags, input identities, gates and complete saved-field errors. The frozen
methods PDFs remain the explicitly labeled **v0.1.0** artifacts. Current function
links are in [FUNCTIONS.md](FUNCTIONS.md); the new control flow is described here.

## Removed work

1. **Discarded initial attempts use an exact stability probe.** A recent initial
   post-chemistry CFL rejection sets a transient hint. On the next first attempt,
   prepare and potential still use the actual half-chemistry state, and the same
   per-face characteristic/volume rates are accumulated in the same order.
   Density reconstruction for discarded transverse fluxes and pressure-source
   work are omitted. The original `min(dt,0.9*CFL/rate)` formula decides whether
   to restore all species/budgets and repeat chemistry at a smaller dt.
2. **A potentially accepted probe is completed.** Full O+ and molecular transport
   runs on that exact prepared state before prediction. Its stability maxima
   must match the probe exactly. A probe is never substituted for an accepted
   flux/source evaluation or reused at a different chemistry state/dt.
3. **Local parallel-face geometry is cached.** Each axis-2 face stores only the
   original weighted mean collision rate, gravity and segment length. These
   values use the original sum/division and are immutable for a mesh epoch.
   Storage covers local face blocks, not a new global canonical allocation.
   It is rebuilt and freed with topology. Ne, pressure and SG coefficients
   remain current-stage quantities.
4. **Post-chemistry response refresh is deferred to prepare.** The two chemistry
   callers return to a driver that unconditionally prepares before the next
   response consumer: RHS/potential, tags/remap or the RHS preceding output.
   Backup, rollback and local reactions use owned molecular state. The old
   immediate halo/response refresh was repeated by prepare with no intervening
   consumer or species update. Every actual consumer still gets current halos
   and conductivity; no physical E-region background is removed or frozen.

`ESF_AUDIT_RATE_PROBE=1` forces initial-attempt probes and also evaluates the
complete transport on the same prepared state. A mismatch of either global
CFL rate or transverse grid-travel rate is fatal. This diagnostic affects work,
not the dt formula or scientific input. Normal mode uses the transient hint;
a restart need not retain it because it selects only which work to perform.

## Same-source-policy benchmark

One verified13350s native/molecular checkpoint pair is advanced to13500s.
Physical roots64x80x64, canonical256x320x64, IGRF-14 and physical L0/L1/L2;
output and checkpoints every150s, including initial/final states. The complete
immediate launch-to-actual-MPI-exit wait includes initialization, cache/restart
loading and output. All runs below reserve256 CPUs on two nodes; both use
64 MPI ranks with four workers each, or256 active threads. No PMPI observer.

| Version in matched R4 allocation | Full wait(s) |
|---|---:|
| Frozen v0.1.0 | 376.018 |
| Current source | 315.240 |

Full elapsed time decreases **16.2%**, with no additional CPUs. This is one
paired launch per version, without a repeated-run uncertainty estimate. It is
not extrapolated to20000s or described as a general strong-scaling result.
The intermediate source-only candidate measured379.702 versus352.799s in its
own separate paired allocation. Those absolute times are history; they are not
subtracted from R4 timings to isolate the refresh edit's contribution.

Both versions accept400 steps with799 attempts and399 initial CFL adjustments;
the reported predictor/corrector failure-retry count remains0. The new source
uses399 stability-probe entries. `rhs_calls=1200` includes them and is not a
count of complete flux evaluations. An accepted probe is completed as described
above. The extra chemistry attempts are retained, not hidden by looser CFL.

The largest saved-field global-amplitude-normalized error is below2e-14.
At both output times the actual physical AMR level maps are byte-identical.
O+ density/momentum, display density, potential and both molecular densities
and momenta are checked, along with true residual and conservation. Momentum
components are normalized independently; this is not pointwise relative error
near vacuum. No equivalence claim is made about unrecorded intermediate maps.

## Verification and scope

All ten configured CTests, eight-rank migration, fresh three-level transport,
fresh two-level night evolution and ASAN/UBSAN four-worker gates pass, with
zero compiler warnings. The fresh night and developed three-level13350-13352s
cases explicitly audit exact stability maxima against full transport. Vendor
MPI/PETSc leak detection is disabled; no ThreadSanitizer claim is made.

The tests reuse a verified immutable checkpoint/cache and preserve every
successful/failed gate and scientific output. Only numerical-source metadata
and integrity checksum in an isolated molecular test copy are rebound after
physical-payload byte checks. Runtime readers remain strict. This harness
conversion is not a general public restart bypass; continue a fresh run using
the same source version, or arrange a separately audited compatible conversion.

The one-dimensional row partition, replicated canonical geometry/molecular
backing and O+ serial face traversal remain important limits. This edit is
qualified for the developed three-level Hong Kong night case, not claimed as
an equal speedup for every dipole, two-dimensional or non-night workflow.
The reduced physical closure and event-calibration limitations of v0.1.0 remain.
Nested phase timers cannot be subtracted or summed into an exclusive budget;
the removed refresh changes the legacy chemistry timer's work scope.

## Root-row layout experiment

A separate same-allocation comparison reuses this exact executable and inputs
to test64MPIx4 (256 active workers) against80MPIx3 (240 workers), with256 CPUs
reserved for both. The physical grid has80 root rows; the latter gives each
rank one physical row. More ranks also replicate more geometry/backing, so
this layout is a measured case-specific option, not an assumed improvement.
The completed [layout evidence](../validation/layout80_qualification.json) records:

| Layout on the same optimized executable/allocation | Active threads | Full wait(s) |
|---|---:|---:|
| 64 MPI x4 | 256 | 309.942 |
| 80 MPI x3 | 240 | 290.586 |

The measured full wait decreases another **6.2%** in this paired allocation.
Both configurations reserve256 CPUs on two nodes. The80-rank affinity probe
verifies disjoint three-CPU domains and three actual workers per rank, and its
short developed-state audit matches complete stability rates exactly. Complete
saved-field errors are below8e-13 and actual AMR maps match byte for byte. The
same400 accepted steps, true residual and conservation checks pass.

For this specific64x80x64 root-grid case,80x3 is the fastest configuration tested
here. It uses more MPI processes/replicated backing but fewer active threads.
Do not apply80 ranks to the smaller teaching grids: rank count must not exceed
the number of physical root rows, and memory/binding must be provisioned.
Intel MPI placement used `-n 80 -ppn 40 -genv I_MPI_PIN_DOMAIN 3 -genv
I_MPI_PIN_ORDER scatter`, `OMP_NUM_THREADS=3`, `OMP_DYNAMIC=false`, and
`OMP_PROC_BIND=false`, on two nodes with256 reserved CPUs. Other MPI vendors
require their own placement settings.

This layout pair uses a separate allocation from the source-change pair.
The source and layout gains are therefore reported separately, not multiplied
or inferred by subtracting absolute timings from different nodes. No production
wall-time extrapolation or asymptotic scaling claim is made.
