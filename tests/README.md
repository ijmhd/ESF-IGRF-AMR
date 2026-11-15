# Test map

Run the configured suite with `ctest --test-dir build --output-on-failure`.
The PETSc/Python configuration has ten tests. Read
[the release evidence](../docs/VALIDATION.md) for which additional workflows
were actually executed on the frozen source; the presence of a script is not
a claim that every optional experiment was rerun for this release.

| CTest name | Source | Main question |
|---|---|---|
| `tube_core` | `test_core.c` | Traced geometry, magnetic/volume and potential kernels |
| `volume_transport` | `test_transport.c` | Conservative mapped-volume transport kernels |
| `electric_driver` | `test_electric_driver.c` | Driver interpolation, coverage and knot clocks |
| `petsc_matrix` | `test_petsc.c` | Distributed tensor algebra versus reference assembly |
| `night_chemistry` | `test_night.c` | Positive local reaction composition, rates and charge ledger |
| `night_multifluid_flux` | `test_night_transport.c` | Molecular mobility/fitted flux and vacuum limits |
| `curved_pdm` | `test_pdm.c` | Physical-volume polynomial moments and bounded reconstruction |
| `collision_background` | `test_collision.c` | Prescribed E population and collision-table responses |
| `spatial_neutral_input` | `test_neutral_grid.c` | Exact spatial binding, dimensions, hash and temperature checks |
| `fixed_amr3_entry_points` | `check_amr.py` | Fixed composite potential, remap and interface transport |

Despite its historical name, `fixed_amr3_entry_points` exercises the two-level
composite backend. The separate fresh three-physical-level gate and developed
three-level restart comparison provide the release's three-level evidence.

The separate MPI migration executable checks arbitrary old/new ownership
intersections, simultaneous cross-peer movement, round trips and exact no-op
state preservation. Launch it explicitly on allocated compute:

```sh
mpiexec -n 8 build/test_amr3_migration
```

Optional retained Python audit drivers:

* `validate.py`: uniform reference-model coupled invariants, MPI comparisons
  and time/geometry refinement experiments. This includes order fitting and
  is not part of the release's performance-only qualification.
* `validate_extensions.py`: uniform reference relaxation/periodic/field-driver
  checks; `--skip-time-refinement` omits its optional order experiment.
* `check_mpi_halos.py`: uniform reference rows/halos across several peers,
  native/PETSc potential paths and periodic/nonperiodic layouts.
* `check_igrf.py`: field evaluator versus an externally supplied exact official
  Fortran source. It separately labels the original and precision-normalized
  references, requires the recorded source hash and a Fortran compiler. No
  reference Fortran solver is bundled or used in ESF evolution.

Use each driver's `--help` and a new work/report path. These optional drivers
compile/execute scientific checks, so use allocated compute. Preserve successful
results and failures. A failed reader after a successful run calls for reading
the saved outputs correctly, not replaying the model to get another reading.
