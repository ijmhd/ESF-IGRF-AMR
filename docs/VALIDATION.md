# Validation and measured performance

Version **v0.1.0** is bound to numerical source ID
`b876245d6051ff517b9121fd6644ecafd762d5ff55ef5c3965895b8edd264c5d`. [Machine-readable evidence](../validation/qualification.json)
records build/executable/input hashes, effective controls, complete saved-field
errors, call counters and completed diagnostics. The tested Linux executable
SHA256 is `3758dca2b4ad8cac095d51ad4b79b34b8bc0a9628fd865d29b0fedb58cbe936c`; a different build prefix/compiler can
produce a different binary while using the same numerical source.

## Completed gates

The allocated Sugon job 29428763 completed with exit0. It built with
GCC12.3.0, Intel MPI2021.3.0, PETSc3.26.1 real-double, OpenMP, Release/O3
and LTO. Compiler warning count: **0**.

* All ten configured CTest tests passed: geometry/volume transport, imposed
  driver, PETSc tensor algebra, local chemistry/flux kernels, curved PDM,
  collision background, strict spatial neutral input and fixed AMR entries.
* An eight-rank ownership-migration test passed with OFI transport.
* Fresh three-level transport and actual night-chemistry workflows passed
  with four workers per rank through0.2s.
* Both fresh workflows passed AddressSanitizer/UndefinedBehaviorSanitizer.
  Vendor-runtime leak detection was disabled (`detect_leaks=0`); this is not
  a leak-sanitizer or ThreadSanitizer qualification.
* The developed-state comparison below passed for all three release variants.

A short fresh workflow checks entry/coupling/I/O, not nonlinear branch growth.
The public CTest suite has ten tests when PETSc and the Python interpreter are
found; the separate migration executable is launched explicitly under MPI.

## Exact timing case

One verified native/molecular pair at13350s is continued to13500s: **150s of
physical advance**, 400 accepted steps, initial/final field outputs and native/
molecular checkpoints every150s. User roots64x80x64, physical L0/L1/L2,
canonical256x320x64, IGRF-14 at2025.0, apex longitude112.25 degrees with
8-degree span, apex250-1600km and conjugate base90km.

The executed case retains PDM4 beta1, CFL0.25, dt maximum2s, potential
request1e-10, the same prescribed background/driver/history, night chemistry,
structure criteria, slab4, buffer1 and quiet60s. Geometry/source comment changes
were bound in **isolated validation copies**, changing only identities, paired
file bindings and integrity checksums. Physical payload bytes were verified
unchanged. Runtime strict readers were not weakened. This is not a public
restart bypass; use a fresh case or a separately audited migration.

The reference is source `54b5c51cc0dc71918e9ffb226aac1fbeefe039f45df19dc870f74bd4480ee124`, the previously
qualified owned-migration/log-cache optimization. Earlier gains against the
older production implementation are not counted again here. Large historical
checkpoints/neutral inputs are not distributed; the documented public examples
run fresh and do not purport to reproduce this large timing case automatically.

## Full launch-to-exit wall time

Every row reserves the same **two-node,256-CPU allocation**. Slurm reports
x86_64,128 CPUs/node and one thread/core. The processor model was not
captured; these absolute times are not a portable hardware benchmark. The active
thread counts differ. Rank domains were checked for disjoint affinity and
sufficient worker placement. Times include initialization, restart/cache
loading, outputs/checkpoints and the actual MPI process exit. No PMPI observer
was preloaded in any row. There is one launch per row, so no repeated-run
uncertainty estimate is available.

| Variant | MPI ranks | Workers/rank | Active threads | Full wait (s) | Reference/time |
|---|---:|---:|---:|---:|---:|
| reference64 | 64 | 1 | 64 | 499.981 | 1.000 |
| optimized64 | 64 | 1 | 64 | 481.027 | 1.039 |
| hybrid32x2 | 32 | 2 | 64 | 537.712 | 0.930 |
| hybrid64x4 | 64 | 4 | 256 | 374.899 | 1.334 |

The same-active-thread64x1 edit reduces measured wall time by
**3.8%**. The fastest measured
configuration is `hybrid64x4` (374.899s), using
256 active threads. Its ratio must not be called a
same-active-core algorithmic speedup. In this grid,32x2 is slower than64x1;
two workers are an example of how to launch, not a blanket recommendation.
The64x4 configuration is useful when enough CPUs are reserved and elapsed time
is the priority. This is not a broad strong-scaling curve or a20000s forecast.

## Numerical comparison and hidden work

At both saved times, all O+ density/momentum, reconstructed display density,
potential, NO+/O2+ density and both molecular momenta are compared with the
reference. Error is `max(abs(q-qref))/max(abs(qref))`, with each momentum
component scaled independently. The largest checked normalized error is
**6.765e-13**, below1e-10. This is global-amplitude normalization, not
pointwise relative error around zero. Actual physical AMR level maps are
byte-identical at both saved frames; no claim is made about unrecorded
intermediate maps. Every completed row independently satisfies O+ mass
error<1e-10 and true potential residual<=1.05e-10.

Release call counters expose400 accepted steps,799 attempts,399 initial
post-chemistry CFL adjustments,1200 RHS calls,1199 chemistry calls and41
tag checks, while reported predictor/corrector failure retries remain0.
Counter scope begins after the initial snapshot and includes the terminal
diagnostic RHS. The step policy is retained; those399 repeats have not been
removed by relaxing a stability check. Nested phase counters cannot be added
into an exclusive100% timeline.

## Remaining limitations and reproducibility

The one-dimensional row partition has only80 physical root rows. Molecular
state and geometry retain full canonical backing; O+ active-cell savings do
not imply proportional memory/speed savings. O+ face/lazy reconstruction and
MPI/PETSc remain on the main thread. These are architectural limits of this
release. Higher active-rank counts or extra workers do not guarantee speedup.

Prescribed temperatures/neutrals, reduced reaction momentum, omitted
photoionization/energy/wind/ion-ion coupling and molecular perpendicular
pressure diffusion limit the science interpretation. Numerical parity does
not calibrate collision multipliers, seed location, event growth time or
Hong Kong/Sanya observations. The methods manuscript states these details.

Rebuild and verify source/text with:

```sh
python3 tools/check_release.py
ctest --test-dir build --output-on-failure
mpiexec -n 8 build/test_amr3_migration
```

Use allocated compute and an MPI rank count allowed by your site. The eight-
rank migration test is small and independent of the science grid. The full
qualification's source, build hashes, exact input identities and scalar
evidence are preserved; no successful scientific run was replayed to repair
a reader or a figure. Production/controller jobs were kept independent.
