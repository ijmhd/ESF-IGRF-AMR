# Geometry and adaptive mesh

## Magnetic and geographic coordinates

The model traces the selected magnetic field in Earth-centred Earth-fixed (ECEF)
metres. IGRF uses radius6371.2km; observation positions use the WGS84 conversion.
Model altitude and apex labels are **spherical** `|x|-6371200m`. A 350km geodetic
HK observation is not obtained by setting spherical height350km at the same
latitude. Map the actual ECEF target using `tube_map`.

The two transverse labels are geographic **apex longitude** and **apex spherical
height**. At a fixed longitude/radius, `tube_apex_label` locates a zero of radial
B. RK4 follows the unit field tangent from the apex to both conjugate base
endpoints. The field coordinate is a weighted trace distance. Its increment is
physical distance divided by an altitude-dependent scale clamped between8km and
80km; it is not uniform arc length, dipole latitude, or geographic latitude.

Each cell is a trilinear hexahedron. Centre metrics supply transverse gradients,
field-aligned geometric tangent and gravity. Two-point tensor Gaussian
quadrature integrates its volume; a separate three-point quadrature provides
cubic reconstruction moments. Folded/zero-volume cells fail. Diagnostics record
alignment, parallel electric metric error and along-tube magnetic-flux jumps.
Replacing dipole by IGRF does not make the coordinates orthogonal; both cross
terms in the conductivity tensor are retained.

## Physical and private levels

| Mode | User physical L0 roots | Internal roots | Canonical geometry |
|---|---|---|---|
| Two levels | Nx x Ny | Nx x Ny | 2Nx x 2Ny |
| Three levels | Nx x Ny | 2Nx x 2Ny | 4Nx x 4Ny |

`ns` stays unchanged. In the three-level implementation private levels `-1,0,1`
mean physical `L0,L1,L2`. A private `-1` macro cell spans2x2 internal roots. The
public output `amr_level` uses physical0/1/2, not these private indices.

`--amr-slab-cells 4` groups tags over four adjacent field segments; it must divide
`ns`. Refinement splits transverse cells by2 in each transverse direction. There
is **no along-field subdivision**, although refinement can vary between slabs.
The full conjugate tube geometry remains present at every physical level.

## Three different notions of mesh

1. **Density DOFs** are actual finite volumes at each k; only their unique owner
   evolves. Nonowner backing slots have zero volume.
2. **Backing/field columns** are the union needed by any refined segment and
   conductivity-driven field request. They may contain many inactive slots.
3. **Canonical cells** provide geometry, physical integration, display sampling
   and the full NO+/O2+ grid. They do not all evolve independent O+ values.

Consequently an electric-field refinement boundary need not follow an O+
depletion contour, and a canonical density image does not prove finest evolution.
Use `amr_level`, activity counts and `structure_frames.csv` for actual density
resolution. At13500s in the audited large case:2219120 active O+ cells,
4739264 backing slots and5242880 canonical molecular cells.

## Structure detector

For each3x3 transverse density stencil, normalize by the larger of local density
and `1e-3 * initial_density_scale`. This is only a detector normalization; it
never changes an evolved density or conductivity.

The raw gradient is `max(|Dx n|,|Dy n|)/norm`; curvature is the Frobenius norm
of the discrete Hessian divided by norm. In `represented` mode, first fit a local
amplitude A of the existing static background-shaped prolongation basis, then
apply the same differences to `(n_neighbor-n_centre)-A*(B_neighbor-B_centre)`.
This does not evolve or subtract an unperturbed control simulation. A bubble
wall/branch remains in the residual; a quiet represented altitude gradient can
avoid unnecessary refinement.

Further criteria:

* Drift curvature: `amr_interval * ||Hessian(u_ExB)|| / local_spacing`.
* Integrated conductivity-tensor curvature requests **field** resolution.
* Prospective-parent information loss compares fine children with the same
  background-preserving conservative prolongation used by remap. The normal
  tolerance is0.005 of the local normalized density.
* Seed marker, if explicitly enabled, protects declared shortest wavelengths.
  The target is32 cells per wavelength and8 cells across the envelope; locations
  unresolved even on the finest mesh are counted. Marker-bearing information
  loss can use the tighter `min(0.005,0.1*seed_amplitude)` tolerance.

Initial and branch wavelengths must be declared when their protection is needed;
AMR cannot reconstruct a perturbation already erased on a coarse initial mesh.
The late-time qualification explicitly uses seed-guard0, because it continues
an existing developed solution rather than inventing a missing historical tracer.

## Buffering and coarsening

Apply a transverse root-cell buffer to structure requests, then promote each
request over its along-field slab. There is no unconditional OR over all k.
The three-level path enforces a refined parent around requested L2 cells,
including neighboring slabs and corner contacts. Field approximation constraints
are built separately from density ownership.

Keep-refined thresholds are lower than refine thresholds. A cell must also stay
quiet and resident for `--amr-quiet-seconds` before coarsening. A tag check occurs
at the configured maximum interval OR when accumulated transverse travel reaches
half an internal root-cell scale. Thus `--amr-interval 30` can produce actual
regrids only3-4s apart in a moving fine-grid case. It is not the observed regrid
period. The qualified policy uses60s quiet time; zero-hold simultaneous coarsening
remains outside the qualified configuration.

## Transfer and shared fluxes

Restriction preserves `sum(V*n)` and `sum(V*m)`. Density prolongation separates
the background-shaped component and bounded residual slopes, with volume-centred
offsets. Velocity prolongation centres its corrections using child particle
weights, so momentum is conserved. Marker transfer follows the same particles.
Remap checks particle, signed momentum and tracer invariants globally.

A coarse/fine contact has a shared physical integrated flux assembled over its
canonical subfaces. Adjacent owner volumes use opposite signs of that same flux.
All levels advance with one accepted dt, so this implementation does not use a
subcycled Berger-Colella flux register. It also does not evolve B or need MHD
constrained transport. Athena++/SWMF informed the conservation/nesting audit;
their source and general block-tree infrastructure are not incorporated.

MPI repartitions whole internal apex-height rows. Molecules move only across
changed old/new ownership intersections; old owned states remain valid through
O+ remapping. Halo widths and minimum row granularity make very high rank counts
inefficient. Geometry and several caches are still replicated.
