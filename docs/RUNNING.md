# Building and running

For the complete Chinese end-to-end workflow, including source download, PETSc
installation, case editing, restart, 3D field-line cuts and movies, see
[从下载源码到模拟和绘图](FROM_SOURCE_TO_MOVIE_ZH.md).

## Dependencies and variants

Use C11, CMake >=3.18, MPI C, pkg-config and real-double PETSc. PETSc must use
the same MPI as `mpicc`/`mpiexec`. OpenMP is optional. Linux/MPI is the qualified
compute environment. A macOS compiler may need an additional OpenMP runtime;
the default build does not request it. Python 3 is for the launcher/checks;
NumPy/matplotlib are optional readers. LaTeX only rebuilds the methods PDF.

```sh
export PKG_CONFIG_PATH=/path/to/petsc/lib/pkgconfig:$PKG_CONFIG_PATH
pkg-config --cflags --libs PETSc
cmake -S . -B build -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_BUILD_TYPE=Release -DESF_MPI=ON -DESF_PETSC=ON
cmake --build build -j 4
ctest --test-dir build --output-on-failure
```

Follow [PETSc's installation documentation](https://petsc.org/release/install/).
For a static PETSc exporting dependencies in `Libs.private`, provide a wrapper
running `pkg-config --static "$@"` through
`-DPKG_CONFIG_EXECUTABLE=/absolute/path/to/wrapper`. Do not discard required
BLAS/LAPACK/MPI/private libraries.

| CMake option | Default | Meaning |
|---|---|---|
| `ESF_MPI` | ON | MPI uniform reference path; needed by AMR/PETSc |
| `ESF_PETSC` | OFF | Constrained AMR and distributed potential solvers |
| `ESF_OPENMP` | OFF | Audited independent kernels; requires MPI/PETSc |
| `ESF_SANITIZE` | OFF | Address/undefined checks propagated through core |
| `CMAKE_BUILD_TYPE` | CMake default | Set Release for performance |
| `CMAKE_INTERPROCEDURAL_OPTIMIZATION` | OFF | Optional supported compiler LTO |

A geometry/reference-only build needs neither MPI nor PETSc:

```sh
cmake -S . -B build-basic -DESF_MPI=OFF -DESF_PETSC=OFF
cmake --build build-basic -j 4
./build-basic/field_eval data/igrf14coeffs.txt 2025 22.3 114.2 350 13
./build-basic/tube_map data/igrf14coeffs.txt 2025 22.3 114.2 350 13 hk-line.csv
```

These utilities report field components/mapping, not bubble evolution. AMR
requires the PETSc build. The static reference solver remains available.

The geometry commands have positional interfaces:

* `field_eval COEFF YEAR LAT_DEG LON_DEG ALT_KM DEGREE` prints north/east/down
  field and its magnitude in nT, evaluated at a WGS84 geodetic position.
* `tube_map COEFF YEAR LAT_DEG LON_DEG ALT_KM DEGREE LINE_CSV` prints a JSON
  target/apex summary and writes the conjugate trace in ECEF, geodetic and
  spherical coordinates. It does not generate an ESF initial state.
* `night_probe O NO O2 R_NO R_O2 ALPHA_NO ALPHA_O2 END DT` integrates only
  local reactions and prints the final species/conversion/recombination
  ledger. The first three values are SI number densities m^-3, the conversion
  rates are s^-1, recombination coefficients m^3/s, and END/DT seconds. These
  are direct coefficients, not temperature/neutral inputs. It has no spatial
  transport, electric potential or photoionization. `night_rates` in `night.c`
  is the separate temperature/neutral-to-coefficient function.

## Explicit cases

The launcher expands checked-in JSON; it prints every effective argument. It
adds only explicit overrides and repository-relative paths. Boolean values are
`0` or `1`; every executable option takes one value, except `--help`.

| Case | Field | Levels | Root grid | Stop | Purpose |
|---|---|---|---|---|---|
| `smoke` | aligned dipole | 3 | 4x4x16 | 0.2 s | Interface/update/I/O |
| `dipole` | aligned dipole | 2 | 16x20x32 | 600 s | Idealized evolution |
| `hk-amr2` | IGRF-14 | 2 | 16x20x32 | 600 s | HK-region geometry |
| `hk-amr3` | IGRF-14 | 3 | 16x20x32 | 600 s | Extra transverse level |
| `hk-night-demo` | IGRF-14 | 2 | 16x20x32 | 600 s | Synthetic-neutral night chemistry |

The night demo starts **2025-03-20 14:00 UTC = 22:00 HKT**, to exercise the
whole-conjugate night guard. The large qualification starts at 12:00 UTC.
Neither clock is 20 MLT. A short coarse run does not establish resolved branches.

```sh
python3 tools/run_case.py hk-amr3 --ranks 4 --out runs/first --dry-run
python3 tools/run_case.py hk-amr3 --ranks 4 --out runs/first
```

For other physics/geometry, edit a copied example argument list and invoke
`esf3d` directly, or deliberately extend the launcher. All controls are in
[PARAMETERS.md](PARAMETERS.md). Use new output directories; existing runs are
preserved. Never infer a physical background from a diagnostic cadence.

### Editing a case explicitly

Each JSON contains a flat `arguments` list of option/value pairs. To keep the
original examples intact, copy one to a separate working directory, edit its
values, and launch it directly. For example, the following *prepares a command*
for a three-level grid with user roots 32x40x48 and a 3600 s stop. It retains
the original case physics and resolves file paths before changing directories:

```python
import json, pathlib, shlex, subprocess
repo = pathlib.Path.cwd().resolve()  # run this from the cloned repository
case = json.loads((repo / "examples/hk-amr3.json").read_text())
args = case["arguments"][:]
for key, value in {"--nx": "32", "--ny": "40", "--ns": "48",
                   "--end": "3600", "--output-every": "150"}.items():
    args[args.index(key) + 1] = value
for key in ["--coeff", "--neutral-profile", "--collision-profile", "--electric-driver"]:
    if key in args:
        i = args.index(key) + 1
        args[i] = str((repo / args[i]).resolve())
out = repo / "runs/custom-hk"
assert not out.exists(), "Choose a new output directory"
args += ["--out", str(out)]
out.parent.mkdir(parents=True, exist_ok=True)
(out.parent / "custom-hk-effective.json").write_text(json.dumps(args, indent=2))
command = ["mpiexec", "-n", "4", str(repo / "build/esf3d"), *args]
print(shlex.join(command))
# After inspecting the printed case, run on allocated compute:
# subprocess.run(command, check=True)
```

Here the canonical grid becomes 128x160x48. A new geometry, neutral-position
mapping and initial state are required; a different grid cannot reuse an old
checkpoint or spatial-neutral file. A height-only teaching neutral table has
different validation from a position-bound spatial table. Increasing `ns`
improves the along-field/E-region sampling separately from transverse AMR.
The example above is a setup recipe, not an additional qualified large run.

Changing `--amr-levels` changes physical levels, while `--amr-mode fixed` holds
an AMR layout and `--amr-mode dynamic` adapts it. `--amr-mode none` uses the
AMR backend without density refinement; omitting the option selects the
separate uniform reference backend. `--amr-mode fine` resolves the available
fine density layout. Consult the grid guide for constraints on each entry.

## MPI/OpenMP and Slurm

MPI owns whole apex-height rows. For three levels, at most user `ny` ranks are
supported. Thin slabs have expensive halos and coarse balancing granularity.
The large qualified three-level process used roughly 4.7-4.9 GiB peak RSS; that
is not an estimate for every grid/PETSc build.

```sh
cmake -S . -B build-omp -DCMAKE_BUILD_TYPE=Release \
  -DESF_MPI=ON -DESF_PETSC=ON -DESF_OPENMP=ON
cmake --build build-omp -j 4
export OMP_NUM_THREADS=2 OMP_DYNAMIC=false
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 tools/run_case.py hk-amr3 --binary build-omp/esf3d \
  --ranks 4 --threads 2 --out runs/hybrid
```

MPI/PETSc calls stay on the main thread; the executable checks thread support.
Setting `OMP_NUM_THREADS` does not adjust CPU binding. Use your MPI/runtime's
placement rules. The Intel MPI qualification uses disjoint rank domains with
at least as many CPUs as workers, and `OMP_PROC_BIND=false` so workers inherit
the domain. Open MPI has different options; do not copy Intel-specific flags.
See [examples/run.slurm](../examples/run.slurm), then set local account, partition
and module paths. Record active threads and reserved CPUs separately.

## Restart

`--end` is the absolute elapsed model stop time, not additional duration.
Night runs require BOTH native and molecular files from a verified same-time
pair. The molecular file binds to the native file hash; mixed epochs fail.

```sh
python3 tools/run_case.py hk-night-demo --ranks 2 --end 1200 \
  --restart runs/night/checkpoint_latest.bin \
  --night-restart runs/night/night_state_latest.bin \
  --out runs/continued
```

Keep physics, grid, seed, AMR thresholds, neutral/driver inputs and source
version consistent. Changing ranks is supported within decomposition limits.
Keep original files. Source edits, including comments in geometry, change strict
fingerprints: use a fresh run or a separately audited conversion, never disable
checks. Geometry caches are native-ABI files, not portable restart states.

## Troubleshooting

* PETSc not found: inspect `PKG_CONFIG_PATH`, MPI consistency and static private
  dependencies. Complex/non-double configurations are unsupported.
* Too many ranks: reduce ranks or deliberately redesign the physical grid.
* Night illumination rejection: check UTC, both conjugate ends and stop time;
  do not bypass the night-only guard.
* Neutral dimensions/positions mismatch: spatial input is bound to exact
  canonical geometry and Ti/Te. Regenerate it, do not reshape its bytes.
* Checkpoint mismatch: inspect pairing, hashes, source and effective parameters.
* Positivity/CFL or true-residual failure: preserve logs and investigate; the
  solver does not silently clip density or accept a failed potential solve.
* `retries=0` is not zero attempted CFL adjustments. Inspect the new
  `performance_counters.json` reason counts.

Use actual frame times and actual AMR labels when plotting. A reconstructed
coarse region has not evolved at finest resolution. Never hold the last saved
field while continuing to advance the displayed movie clock.

## Plotting a complete saved slice

Install the optional reader packages from `requirements-plot.txt`. On allocated
compute holding the saved fields:

```sh
python3 tools/plot_slice.py runs/first --out figures/first-slice
python3 tools/plot_slice.py runs/first --time 600 --field-index 16 \
  --quantity electron --out figures/first-t600
```

The utility chooses the latest complete frame unless `--time` is specified and
writes PDF, PNG and input-hash provenance. Use a new output prefix. Its axes are
apex longitude and apex altitude at a fixed field-label index, **not** geographic
350 km or an equatorial plane. White outlines mark changes of actual physical
density refinement level; the second panel labels L0/L1/L2 explicitly. Full
geographic/ISR sampling must use ECEF geometry, as described in the output guide.
