# Output and restart contracts

## Frames and units

Read `metadata.json`, `frames.csv` AND `structure_frames.csv` first. Each frame
binds to its mesh epoch; do not decode a later dynamic frame with initial
`leaves.csv`. Canonical indices are [j,i,k], k contiguous. Use binary files on a
compatible IEEE-double/little-endian platform or decode explicitly. Portable
checkpoint encoding is distinct from native visualization/cache files.

| File | Shape / meaning |
|---|---|
| `state_NNNN.bin` | float64[2,nleaf,ns]: O+ normalized density and N*u_parallel |
| `leaves_EEEE.csv` | id, private backing level, i,j, internal root |
| `leaf_metrics_EEEE.bin` | float64[nleaf,ns,6]: volume,bg,nu,gravity,area_lo,area_hi; use metadata order |
| `display_density_NNNN.bin` | float64[fine_ny,fine_nx,ns]: conservative reconstructed O+ only |
| `phi_NNNN.bin` | float64[fine_ny+1,fine_nx+1]: canonical potential in volts |
| `molecular_NNNN.bin` | float64[2,fine_ny,fine_nx,ns]: NO+,O2+ normalized density |
| `molecular_momentum_NNNN.bin` | same shape: each normalized density times own parallel velocity |
| `amr_level_NNNN.bin` | uint8[fine_ny,fine_nx,ns]: actual physical density level0/1/2 |
| `density_mask_EEEE.bin` | uint8[internal_ny,internal_nx,ns]: actual segmented L1/L2 requests |
| `tag_reasons_NNNN.bin` | native uint16[internal_ny,internal_nx,ns]: reason bit mask |
| `tag_values_NNNN.bin` | float64[internal_ny,internal_nx,ns,6]: G,C,V,loss,rawG,rawC |
| `marker_NNNN.bin` | optional normalized material quantity on backing density state |
| `nodes.bin` | ECEF geometry; consult metadata for exact node shape |

Multiply stored normalized densities by1e12 to get m^-3. Where density is
positive, divide momentum by density for m/s; exactly zero density requires
zero momentum. `display_density` is O+, not electron density. With night mode,
canonical Ne is reconstructed O+ plus NO+ plus O2+. At coarse O+ cells, displayed
canonical variation is reconstruction, not an independently evolved solution.

The initial legacy `leaf_metrics.bin` writes volume,bg,nu,gravity,area_lo,area_hi;
epoch-specific metrics use that same explicitly written order. Internal A3Metric
member order is different: never use a raw C struct cast as the file schema.

## Scalar diagnostics

`diagnostics.csv` records mass balance corrected for integrated boundary flux
and O+ conversion loss, actual potential true residual, CFL rate and cumulative
potential calls. `night_diagnostics.csv` records conversion/recombination,
molecular mass, charge error and total-ion reaction/boundary balance.
These mass/particle ledgers integrate the STORED normalized density: multiply
absolute particle totals or time-integrated particle fluxes by1e12 for physical
particle numbers. They are not kilograms. Relative balances are independent of
this factor; physical O+ momentum additionally needs its ion mass.
`regrid.csv` records actual changed layouts, ownership movement, transfer errors
and potential change. A changed approximation space need not leave phi identical.

`completion.json` summarizes the completed interval. Its field `cells` is backing
slots, not actual activity. Obtain activity from `three_level_status.json` /
`structure_status.json`. Completion phase timers are rank-local and overlap;
regrid includes nested prepare/solve. `chemistry_wall_s` also includes molecular
response refresh/halo waiting. Do not add these into an exclusive time budget.
Full synchronous launch-to-exit timing is the end-to-end performance measure.

`performance_counters.json` starts after the initial snapshot and includes the
terminal diagnostic RHS. It records RHS/tag/chemistry calls, total attempted
steps, initial post-chemistry CFL adjustments and reported failure retries. The
last are distinct concepts. Historical step/solve counts may include restart
history. `reported_failure_retries` and `segment_failure_retries` exclude it;
`start_retries` records the restored count and `retries_total` retains the full
history, also retained by `completion.json.retries` and the next checkpoint.

## Native and molecular checkpoints

New native AMR checkpoints begin with `ESFSRC01` (8 bytes), then the 64 ASCII
hex characters of the complete compiled numerical-source ID (`ESF_NIGHT_ID`,
the hash of all src C/header and include header hashes). The existing native
payload follows; its checksum now also covers this source envelope. This
binding applies to non-night as well as night runs. A mismatching source ID is
fatal before the state is accepted. The original unbound checkpoint formats,
including historical static-state imports, must be continued with their
original executable; this version does not silently migrate or rebind them.
Keep the original version, inputs and results. No bypass flag is provided.

The inner native magic distinguishes ESFAMR01 legacy, ESFAMR02 structure,
ESFAMR03 driver extension, ESFAMR04 spatial neutrals, and ESFAMR05 three physical
levels. All have checksummed little-endian words/doubles. Structure formats
contain geometry source/input signatures, physical/numerical parameter arrays,
union/density/field masks, quiet/birth clocks, optional marker, density/momentum,
potential and original global conservation histories.

The ESFNIT02 file has an8-byte magic; seven uint64 values describe dimensions,
cell count, neutral input hash, paired native-file hash and algorithm version;
a64-byte numerical-source ID follows. Nine doubles store time/UTC/temperatures/
reaction and boundary histories, then the two ion-density and two momentum
arrays, and an8-byte FNV checksum. The pair is complete only when both atomically
renamed files match. A run can stop between the two writes, so verify binding
before treating latest filenames as an archived restart.

FNV is an integrity check, not a cryptographic security boundary. Release
manifests additionally record SHA256. Runtime readers reject mismatch/nonfinite
data; they do not silently reinterpret a checkpoint.

## Geometry cache and spatial neutral input

ESFGEO01 caches native-ABI geometry with shape/parameter/input/code/ABI checks and
a payload checksum. It is not portable across ABI/source changes and is not a
restart state. Spatial ESFNEU01 is a strict little-endian native-double input:
8-byte magic, seven uint64 fields (version,nx,ny,ns,6,position_hash,value_hash),
Ti/Te doubles, then float64[ny,nx,ns,6] containing O,N2,O2 number densities m^-3,
Tn K, O+ collision s^-1 and electron collision s^-1. Positions must hash exactly
to canonical cell centres. A height neutral table instead has columns
height_m,O_m3,N2_m3,O2_m3,Tn_K and positive log-interpolated values.

Reader implementation and every field/order are available in the annotated
checkpoint/neutral modules and the source atlas. Avoid treating visualization
binaries, geometry cache, neutral input and portable checkpoint as one format.
